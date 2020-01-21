###
# PRISM REPURPOSING SECONDARY DATASET DATA PROCESSING SCRIPT
# Author : Mustafa Anil Kocak, mkocak@broadinstitute.org
# Last Modified: Nov 21, 2019
# 
# For the raw and processed data please visit https://www.depmap.org
# For the detailed description of the dataset : https://doi.org/10.1101/730119
#
# The script below includes all the data processing steps for the secondary screen dataset.
###

# ------
# LOAD THE REQUIRED LIBRARIES
# ----

## Please load the missing libraries by uncommenting the first three lines below
#install.packages("tidyverse", "data.table", "magrittr", "reshape2", "drc")
#install.packages("BiocManager")
#BiocManager::install("sva")

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(sva)
library(drc)

# ----
# AUXILARY FUNCTIONS 
# ----

compute_auc = function(l, u, ec50, h, md, MD) {
  #if( l > 1) l = 1
  #if( l < 0) l = 0
  f1 = function(x) pmax(pmin((l + (u - l)/(1 + (2^x/ec50)^h)),1),0)
  integrate(f1, log2(md),log2(MD))$value/(log2(MD/md))
}

compute_log.ic50 = function(l, u, ec50, h, md, MD) {
  if((l >= 0.5) | (u <= 0.5)){
    return(NA)
  }else{
    f1 = function(x) (l + (u - l)/(1 + (2^x/ec50)^h) - 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root, error = function(x) NA))
  } 
}


apply_combat <- function(Y) {
  df <- reshape2::dcast(Y, row_name + screen_id + pool_id ~ detection_plate, value.var = 'log2.Viability', fun.aggregate = function(x) median(x, na.rm = T)) %>%
    #tidyr::unite(cond, screen_id, pool_id) 
    tidyr::gather(key = "detection_plate", value = "value", 4:dim(.)[2]) %>%
    drop_na() %>%
    tidyr::unite(cond, screen_id, pool_id, detection_plate, remove = F) 
  
  batch <- df$cond
  m <- rbind(df$value, rnorm(length(df$value), mean =  mean(df$value, na.rm = T), sd = sd(df$value, na.rm = T)))
  
  combat <- tryCatch(sva::ComBat(dat = m, batch = batch) %>%
                       t() %>%
                       as.data.frame() %>%
                       dplyr::mutate(row_name = df$row_name, detection_plate = df$detection_plate, 
                                     screen_id = df$screen_id, 
                                     pool_id = df$pool_id) %>%
                       dplyr::rename(log2.Viability.cb = V1) %>%
                       dplyr::select(-V2), 
                     error = function(e) dplyr::rename(df, log2.Viability.cb = value))
  
  
  Y %>% 
    dplyr::left_join(combat) %>%
    .$log2.Viability.cb
}

# -----
# SPECIFY THE FILE PATHS
# ----

## Please modify the following three lines to the appropriate adresses to the files downloaded from depmap.org
treatment_info_path = "treatment_info.csv"
pooling_info_path = "pooling_info.csv"
secondary_MFI_path = "MFI.csv"

# ----
# LOAD THE DATA
# ----

treatment_info = data.table::fread(treatment_info_path) %>%
  dplyr::select(column_name, perturbation_type, screen_id, detection_plate, compound_plate, well, broad_id, dose)

pooling_info = data.table::fread(pooling_info_path) %>%
  dplyr::select(row_name, screen_id, detection_pool, pool_id)

DATA = data.table::fread(secondary_MFI_path) %>%
  tidyr::gather(key = "column_name", value = "MFI", -V1) %>%
  dplyr::rename(row_name = V1) %>% 
  dplyr::left_join(treatment_info) %>%
  dplyr::left_join(pooling_info) %>%
  dplyr::filter(is.finite(MFI)) %>%
  dplyr::mutate(log2.MFI = log2(MFI))

# ----
# COMPUTE THE OUTLIER POOLS AND SSMD'S FOR EACH CELL LINE
# ----

ARTIFACTS = DATA %>%
  dplyr::group_by(screen_id, row_name, detection_plate) %>%
  dplyr::mutate(log2.MFI = log2.MFI - median(log2.MFI, na.rm = T)) %>% 
  dplyr::group_by(screen_id, detection_plate, well, pool_id, perturbation_type, compound_plate) %>%
  dplyr::summarize(median.lfc = median(log2.MFI, na.rm = T)) %>%
  dplyr::group_by(screen_id, well, compound_plate) %>%
  dplyr::mutate(diff.median.lfc = (median.lfc - median(median.lfc, na.rm=T)) / mad(median.lfc, na.rm=T)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(is_outlier = abs(diff.median.lfc) > 5)

SSMD_TABLE = DATA %>% 
  dplyr::left_join(ARTIFACTS) %>% 
  dplyr::filter(!is_outlier) %>%
  dplyr::select(-median.lfc, -diff.median.lfc, -is_outlier) %>%  
  dplyr::group_by(screen_id, detection_plate, compound_plate, row_name) %>%
  dplyr::summarise(neg.con.med = median(log2.MFI[perturbation_type == "vehicle_control"]),
                   neg.con.mad = mad(log2.MFI[perturbation_type == "vehicle_control"]),
                   pos.con.med = median(log2.MFI[perturbation_type == "positive_control"]),
                   pos.con.mad = mad(log2.MFI[perturbation_type == "positive_control"])) %>% 
  dplyr::mutate(ssmd = (neg.con.med - pos.con.med) / sqrt(neg.con.mad^2 + pos.con.mad^2))

# ----
# FILTER AND PROCESS THE DATA
# ----

DATA = DATA %>% 
  dplyr::left_join(SSMD_TABLE) %>%
  dplyr::filter((is.finite(ssmd) & (ssmd >= 2)) | (pool_id == "CP01")) %>%
  dplyr::left_join(ARTIFACTS) %>%
  dplyr::filter((abs(diff.median.lfc) <= 5) | (pool_id == "CP01"))%>% 
  dplyr::select(-ssmd, -neg.con.med, -neg.con.mad, -pos.con.med, -pos.con.mad, -median.lfc, -diff.median.lfc) %>%
  dplyr::group_by(column_name) %>%
  dplyr::mutate(log2.nMFI = log2.MFI - median(log2.MFI[pool_id == "CP01"])) %>%
  dplyr::mutate(log2.nMFI = ifelse(is.finite(log2.nMFI), log2.MFI, NA)) %>% 
  dplyr::group_by(row_name, detection_plate, screen_id) %>%
  dplyr::mutate(log2.Viability = log2.nMFI - median(log2.nMFI[perturbation_type == "vehicle_control"], na.rm = T)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(pool_id != "CP01") %>% 
  tidyr::unite(condition, screen_id, compound_plate, broad_id, dose, remove = F) %>% 
  tidyr::unite(condition2, screen_id, compound_plate, well, remove = F) %>% 
  dplyr::mutate(condition = ifelse(perturbation_type == "experimental_treatment", condition, condition2)) %>%
  dplyr::select(-condition2) %>%
  split(.$condition) %>%
  purrr::map_dfr(~dplyr::mutate(.x, log2.Viability.cb = apply_combat(.))) 

#----
# MEDIAN-COLLAPSE THE REPLICATES TO GET A SINGLE PROFILE
#----

DATA.median.collapsed = DATA %>% 
  dplyr::filter(perturbation_type == "experimental_treatment") %>%
  dplyr::group_by(row_name, perturbation_type, screen_id, compound_plate, broad_id, dose) %>%
  dplyr::summarise(log2.Viability = median(log2.Viability, na.rm = T), 
                   log2.Viability.cb = median(log2.Viability.cb, na.rm = T)) %>%
  tidyr::unite(col_name, broad_id, dose, screen_id, compound_plate, sep = "::")



# ----
# COMPUTE THE DOSE-RESPONSE PARAMETERS
# ----


DRC_TABLE = DATA %>%
  dplyr::filter(perturbation_type == "experimental_treatment", pool_id != "CP01") %>% 
  dplyr::distinct(row_name, broad_id, dose) %>%
  dplyr::count(row_name, broad_id) %>%
  dplyr::filter(n > 4) %>%
  # sample_n(100) %>%
  dplyr::mutate(ix = 1:n()) 


DRC = tibble()

for(jx in 1:nrow(DRC_TABLE)){
  print(paste0(jx, " of ", nrow(DRC_TABLE)))
  data = DRC_TABLE %>%
    dplyr::filter(ix == jx) %>% 
    dplyr::left_join(DATA) %>%
    dplyr::filter(is.finite(log2.Viability.cb))
  
  fit <- tryCatch(drc::drm(2^log2.Viability.cb ~ dose,
                           data=data,
                           fct=drc::LL.4(fixed=c(NA, NA, 1, NA)),
                           robust="median", na.action = na.omit), 
                  error = function(e) NULL)
  
  #plot(fit, type = "all")
  
  if(!is.null(fit)){
    if(fit$fit$convergence){
      pred <- predict(fit, newdata=data.frame(fit$dataList$dose))
      true <- fit$dataList$origResp
      
      temp <- dplyr::tibble("ix" = jx, 
                            "UpperLimit" = 1, 
                            "EC50" = fit$parmMat[3], 
                            "Slope" = fit$parmMat[1],
                            "LowerLimit" = fit$parmMat[2]) %>%
        dplyr::mutate(auc = compute_auc(LowerLimit, UpperLimit, 
                                        EC50, Slope, 
                                        min(data$dose), max(data$dose)),
                      log2.ic50 = compute_log.ic50(LowerLimit, UpperLimit, 
                                                   EC50, Slope, 
                                                   min(data$dose), max(data$dose)),
                      R2 = 1 - (sum((true - pred)^2, na.rm = T) / sum((true - mean(true))^2, na.rm = T)))
      
      DRC %<>% 
        bind_rows(temp)
    }
  }
}


DRC_TABLE = DRC %>%
  dplyr::left_join(DRC_TABLE) %>% 
  dplyr::select(-ix, -n) 

# ----
# SAVE THE PROCESSED FILES 
# ----

ARTIFACTS %>%
  write_csv("secondary_outlier_pools.csv")

SSMD_TABLE %>%
  write_csv("secondary_ssmd_table.csv")

DATA %>% 
  reshape2::acast(row_name ~ column_name, value.var = "log2.Viability.cb") %>%
  write.csv("secondary_logfold_change.csv")

DATA.median.collapsed %>%
  reshape2::acast(row_name ~ col_name, value.var = "log2.Viability.cb", fun.aggregate = function(x) mean(x, na.rm = T)) %>%
  write.csv("secondary_replicate_collapsed_logfold_change.csv")

DRC_TABLE %>%
  write_csv("secondary_drc_table.csv")
