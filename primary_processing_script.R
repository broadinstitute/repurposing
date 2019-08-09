###
# PRISM REPURPOSING PRIMARY DATASET DATA PROCESSING SCRIPT
# Author : Mustafa Anil Kocak, mkocak@broadinstitute.org
# Last Modified: Aug 9, 2019
# 
# For the raw and processed data please visit https://www.depmap.org
# For the detailed description of the dataset : https://doi.org/10.1101/730119
#
# The script below includes all the data processing steps for the primary screen dataset.
###

# ------
# LOAD THE REQUIRED LIBRARIES
# ----

## Please load the missing libraries by uncommenting the first three lines below
#install.packages("tidyverse", "data.table", "magrittr", "reshape2")
#install.packages("BiocManager")
#BiocManager::install("sva")

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(sva)

# ----
# AUXILARY FUNCTION TO APPLY COMBAT
# ----

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
treatment_info_path = "primary_replicate_treatment_info.csv"
pooling_info_path = "pooling_info.csv"
primary_MFI_path = "primary_MFI.csv"

# ----
# LOAD THE DATA
# ----

treatment_info = data.table::fread(treatment_info_path) %>%
  dplyr::select(column_name, perturbation_type, screen_id, detection_plate, compound_plate, well, broad_id, dose)

pooling_info = data.table::fread(pooling_info_path) %>%
  dplyr::select(row_name, screen_id, detection_pool, pool_id)

DATA = data.table::fread(primary_MFI_path) %>%
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
  dplyr::group_by(screen_id, compound_plate, row_name) %>%
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
  dplyr::filter(is.finite(ssmd), ssmd >= 2) %>%
  dplyr::left_join(ARTIFACTS) %>%
  dplyr::filter(abs(diff.median.lfc) <= 5) %>% 
  dplyr::select(-ssmd, -neg.con.med, -neg.con.mad, -pos.con.med, -pos.con.mad, -median.lfc, -diff.median.lfc) %>% 
  dplyr::group_by(row_name, detection_plate, screen_id) %>%
  dplyr::mutate(log2.Viability = log2.MFI - median(log2.MFI[perturbation_type == "vehicle_control"], na.rm = T)) %>%
  dplyr::ungroup() %>%
  tidyr::unite(condition, screen_id, compound_plate, well, remove = F) %>% 
  split(.$condition) %>%
  purrr::map_dfr(~dplyr::mutate(.x, log2.Viability.cb = apply_combat(.))) 

#----
# MEDIAN-COLLAPSE THE REPLICATES TO GET A SINGLE PROFILE
#----

DATA.median.collapsed = DATA %>% 
  dplyr::filter(perturbation_type == "experimental_treatment") %>%
  dplyr::group_by(row_name, perturbation_type, screen_id, compound_plate, well, broad_id, dose) %>%
  dplyr::summarise(log2.Viability = median(log2.Viability, na.rm = T), 
                   log2.Viability.cb = median(log2.Viability.cb, na.rm = T)) %>%
  tidyr::unite(col_name, broad_id, dose, screen_id, sep = "::")

# ----
# SAVE THE PROCESSED FILES 
# ----

ARTIFACTS %>%
  write_csv("primary_outlier_pools.csv")

SSMD_TABLE %>%
  write_csv("primary_ssmd_table.csv")

DATA %>% 
  reshape2::acast(row_name ~ column_name, value.var = "log2.Viability.cb") %>%
  write.csv("primary_logfold_change.csv")

DATA.median.collapsed %>%
  reshape2::acast(row_name ~ col_name, value.var = "log2.Viability.cb") %>%
  write.csv("primary_replicate_collapsed_logfold_change.csv")
