library(tidyverse)
library(parallel)
library(dr4pl)
library(drc)


# Set the path for the files downloaded from depmap.org/repurposing
path <- "~/Downloads/secondary portal update/"

# Load the data 
cell_line_info <- data.table::fread(paste0(path, "secondary-screen-cell-line-info.csv"))
inst_info <- data.table::fread(paste0(path, "secondary-screen-replicate-treatment-info.csv"))
log2_fold_change <- data.table::fread(paste0(path, "secondary-screen-logfold-change.csv")) %>% 
  column_to_rownames("V1") %>% 
  as.matrix()

  
# Formatting viability matrix and condition annotations
ViabilityConditions <- inst_info %>% 
  dplyr::filter(perturbation_type == "experimental_treatment",
                column_name %in% colnames(log2_fold_change)) %>% 
  dplyr::rename(Label = column_name,
                SampleID = broad_id,
                Dose = dose,
                CompoundPlate = compound_plate) %>%
  dplyr::mutate(DoseUnit = "uM") %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit, CompoundPlate, detection_plate) %>% 
  dplyr::group_by(SampleID, Dose, CompoundPlate) %>%  
  dplyr::arrange(detection_plate) %>% 
  dplyr::mutate(Replicate = 1:n()) %>% 
  dplyr::select(-detection_plate) %>% 
  dplyr::ungroup()
  
ViabilityMatrix <- 2^log2_fold_change[is.na(word(rownames(log2_fold_change), 2, sep = fixed("_"))), unique(ViabilityConditions$Label)]

collapsed_lfc <- ViabilityMatrix %>% 
  reshape2::melt() %>%
  dplyr::filter(is.finite(value)) %>% 
  dplyr::rename(Label = Var2, FC = value) %>% 
  dplyr::left_join(ViabilityConditions) %>% 
  dplyr::mutate(LFC = log2(FC)) %>% 
  dplyr::group_by(SampleID, Dose, DoseUnit, CompoundPlate, Var1) %>% 
  dplyr::filter(n() > 1) %>%
  dplyr::summarise(LFC = median(LFC)) %>% 
  dplyr::ungroup()

selected_profiles <- collapsed_lfc %>% 
  dplyr::count(SampleID, CompoundPlate) %>% 
  dplyr::group_by(SampleID) %>% 
  dplyr::top_n(1, n) %>%
  dplyr::ungroup()

collapsed_lfc <- selected_profiles %>% 
  dplyr::left_join(collapsed_lfc) %>% 
  dplyr::left_join(inst_info %>% 
                     dplyr::distinct(broad_id, name) %>%
                     dplyr::rename(SampleID = broad_id)) %>% 
  dplyr::mutate(Label = paste0(toupper(name), "(", SampleID, ") @",Dose, " ", DoseUnit)) %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit, Var1, LFC)

ViabilityCollapsedConditions <- collapsed_lfc %>% 
  dplyr::distinct(Label, SampleID, Dose, DoseUnit)

ViabilityCollapsedMatrix <- 2^reshape2::acast(collapsed_lfc, Var1 ~ Label, value.var = "LFC")

ViabilityConditions <- ViabilityConditions %>% 
  dplyr::semi_join(selected_profiles) 

ViabilityMatrix <- ViabilityMatrix[, unique(ViabilityConditions$Label)]

# Fitting Dose Response Curves-----

compute_auc <- function(LL, UL, Inflection, Slope, md, MD) {
  f1 = function(x) pmax(pmin((UL + (LL - UL)/(1 + (2^x/Inflection)^Slope)), 1, na.rm = T), 0, na.rm = T)
  return(tryCatch(integrate(f1, log2(md), log2(MD))$value/(log2(MD/md)),
                  error = function(e) {print(e); NA}))
}
compute_log_ic50 <- function(LL, UL, Inflection, Slope, md, MD) {
  if((LL >= 0.5) | (UL <= 0.5)) {
    return(NA)
  } else {
    f1 = function(x) (UL + (LL - UL)/(1 + (2^x/Inflection)^Slope)- 0.5)
    return(tryCatch(uniroot(f1, c(log2(md), log2(MD)))$root,
                    error = function(x) NA))
  }
}


LFC.parallel <- ViabilityMatrix %>% 
  log2() %>% 
  reshape2::melt() %>% 
  dplyr::rename(depmap_id = Var1,
                Label = Var2,
                LFC_cb = value) %>% 
  dplyr::filter(is.finite(LFC_cb)) %>% 
  dplyr::left_join(ViabilityConditions) %>% 
  tidyr::unite(px, SampleID, CompoundPlate, remove = FALSE) 


LFC.parallel <- lapply(unique(LFC.parallel$px), function(x) dplyr::filter(LFC.parallel, px == x))

# Encapsulating get_best_fit and compute_MSE_MAD functions for parallellization
f <- function(df) { 
  require(tidyverse)
  compute_MSE_MAD <- function(FC, dose,  UL, LL,  Slope, Inflection) {
    FC.pred = UL  + (LL -UL )/(1 + (dose/Inflection)^Slope)
    residuals = FC - FC.pred
    return(list(mse = mean(residuals^2), mad = median(abs(residuals))))
  }
  get_best_fit <- function(FC, dose, UL_low=0.8, UL_up=1.01, slope_decreasing=TRUE) {
    require(dr4pl)
    require(drc)
    require(tidyverse)
    require(magrittr)
    
    # Fits a number of alternate models  to the DRC and chooses the best fit.
    
    # UL low is the lowerbound of UL we pass to the optimizer and UL_up is the upper bound of UL that we pass to the optimizer
    # fomat of output will be:-
    # results.df <- data.frame("fit_name"=character(),"Lower_Limit"=double(),
    #                          "Upper_Limit"=double(), 
    #                          "Slope"=double(),
    #                          "Inflection"=double(), 
    #                          "MSE"=double(), "MAD" =double(),
    #                          "frac_var_explained"=double())
    
    
    
    riemann_AUC <- mean(pmin(1,FC)) ## mean fold-change after rounding FC to 1.
    var_data = var(FC)
    slope_bound <- ifelse(slope_decreasing, 1e-5, Inf)  # bound the slopes by default unless passed another option
    
    
    results.df <- list(); ix = 1
    
    
    # FIT 1 ---------
    drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                    fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                    lowerl = c(-slope_bound,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)), # ?? 
                           error = function(e)
                           {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
    # "slope" in drc package is -ve of slope in dr4pl package
    
    
    if (drc_model$fit$convergence){
      mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
      
      results.df[[ix]] <- tibble( fit_name = "drc_drm_constrained",
                                  Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                  Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                  Slope = -as.numeric(drc_model$coefficients[[1]]),
                                  Inflection = as.numeric(drc_model$coefficients[[4]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # FIT 2 ---------
    drc_model <-  tryCatch(drc::drm(FC ~ dose, data= data.frame(FC = FC, dose = dose),
                                    fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")),
                                    lowerl = c(-Inf,0.0, UL_low, -Inf),upperl = c(Inf,1.01,UL_up, Inf)), # ?? 
                           error = function(e)
                           {return(list(convergence=FALSE, error=TRUE,fit=list(convergence=FALSE)))})
    # "slope" in drc package is -ve of slope in dr4pl package
    
    
    if (drc_model$fit$convergence){
      mse_mad <- compute_MSE_MAD(FC, dose, as.numeric(drc_model$coefficients[[3]]), as.numeric(drc_model$coefficients[[2]]),
                                 -as.numeric(drc_model$coefficients[[1]]), as.numeric(drc_model$coefficients[[4]]))
      # "slope" in drc package is -ve of slope in dr4pl package and so -ve sign needs to be put in here.
      
      results.df[[ix]] <- tibble( fit_name = "drc_drm_unconstrained",
                                  Lower_Limit = as.numeric(drc_model$coefficients[[2]]),
                                  Upper_Limit = as.numeric(drc_model$coefficients[[3]]),
                                  Slope = -as.numeric(drc_model$coefficients[[1]]),
                                  Inflection = as.numeric(drc_model$coefficients[[4]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # FIT 3 ------
    dr4pl_initMan_optNM <- tryCatch(dr4pl(dose, FC,
                                          init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_2 = 8*min(dose),
                                                                         theta_3= -3, theta_4 = 0.01),
                                          lowerl = c(UL_low, -Inf, -Inf, 0),
                                          upperl = c(UL_up, Inf, slope_bound, 1.01),
                                          method.optim="Nelder-Mead"),
                                    error= function(e){return(list(convergence=FALSE, error=TRUE))}
    )
    
    if (dr4pl_initMan_optNM$convergence==FALSE){
      if (!is.null(dr4pl_initMan_optNM$dr4pl.robust)) {
        dr4pl_initMan_optNM <- dr4pl_initMan_optNM$dr4pl.robust
      }
    }
    
    if (dr4pl_initMan_optNM$convergence){
      mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_initMan_optNM$parameters[[1]], dr4pl_initMan_optNM$parameters[[4]],
                                 dr4pl_initMan_optNM$parameters[[3]], dr4pl_initMan_optNM$parameters[[2]])
      
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initMan_constrained_optNM",
                                  Lower_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[4]]),
                                  Upper_Limit = as.numeric(dr4pl_initMan_optNM$parameters[[1]]),
                                  Slope = as.numeric(dr4pl_initMan_optNM$parameters[[3]]),
                                  Inflection = as.numeric(dr4pl_initMan_optNM$parameters[[2]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # FIT 4 -----
    dr4pl_unconstrained <- tryCatch(dr4pl(dose, FC,
                                          init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.3),
                                          method.init = "logistic",
                                          lowerl = c(0.99, -Inf, -Inf, 0),
                                          upperl = c(1.01, Inf, Inf, 1.01)),
                                    error = function(e) {print(e); return(NA)})
    
    if (!all(is.na(dr4pl_unconstrained))) {
      if (!dr4pl_unconstrained$convergence) {
        dr4pl_unconstrained <- dr4pl_unconstrained$dr4pl.robust
      }
    }
    
    
    param <- tryCatch(dr4pl_unconstrained$parameters, error = function(e) return(NA))
    if (!all(is.na(param))){
      if(as.numeric(dr4pl_unconstrained$parameters[[3]])<slope_bound){ ### while slope bound is not passed to this last optimizer, we do not accept a solution not within the bound
        mse_mad <- compute_MSE_MAD(FC, dose, dr4pl_unconstrained$parameters[[1]], dr4pl_unconstrained$parameters[[4]],
                                   dr4pl_unconstrained$parameters[[3]], dr4pl_unconstrained$parameters[[2]])
        results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_unconstrained",
                                    Lower_Limit = as.numeric(dr4pl_unconstrained$parameters[[4]]),
                                    Upper_Limit = as.numeric(dr4pl_unconstrained$parameters[[1]]),
                                    Slope = as.numeric(dr4pl_unconstrained$parameters[[3]]),
                                    Inflection = as.numeric(dr4pl_unconstrained$parameters[[2]]),
                                    MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
        ix = ix + 1
      }
    }
    
    # FIT 5 ----
    dr4pl_initL <- tryCatch(dr4pl(dose, FC,
                                  init.parm = dr4pl::dr4pl_theta(theta_1 = 1, theta_4 = 0.005),
                                  method.init = "logistic",
                                  lowerl = c(UL_low, -Inf, -Inf, 0),
                                  upperl = c(UL_up, Inf, slope_bound, 1.01)),
                            error= function(e){return(list(convergence=FALSE, error=TRUE))}
    )
    
    if (dr4pl_initL$convergence==FALSE){
      if (!is.null(dr4pl_initL$dr4pl.robust)) {
        dr4pl_initL <- dr4pl_initL$dr4pl.robust
      }
    }
    
    if (dr4pl_initL$convergence){
      mse_mad <- compute_MSE_MAD(FC,dose, dr4pl_initL$parameters[[1]], dr4pl_initL$parameters[[4]],
                                 dr4pl_initL$parameters[[3]], dr4pl_initL$parameters[[2]])
      
      results.df[[ix]] <- tibble( fit_name = "dr4pl_initL_constrained",
                                  Lower_Limit = as.numeric(dr4pl_initL$parameters[[4]]),
                                  Upper_Limit = as.numeric(dr4pl_initL$parameters[[1]]),
                                  Slope = as.numeric(dr4pl_initL$parameters[[3]]),
                                  Inflection = as.numeric(dr4pl_initL$parameters[[2]]),
                                  MSE = mse_mad$mse, MAD = mse_mad$mad, frac_var_explained = 1-mse_mad$mse/var_data)
      ix = ix + 1
    }
    
    # Choose the best fit among the successful fits------
    results.df <- dplyr::bind_rows(results.df) 
    
    # if (nrow(results.df)>0){
    #   results.df <-  dplyr::filter(results.df, frac_var_explained > 0)
    # }
    
    if (nrow(results.df)>0){
      results.df <- results.df %>%
        dplyr::arrange(desc(frac_var_explained)) %>% 
        head(1) %>% 
        dplyr::mutate(successful_fit = TRUE, AUC_Riemann = as.numeric(riemann_AUC) ) 
    }else{
      results.df  <- data.frame(successful_fit=FALSE, AUC_Riemann = riemann_AUC) 
    }
    
    return (results.df)
  }
  
  df %>%
    dplyr::group_by(SampleID, CompoundPlate, depmap_id) %>%
    dplyr::summarise(get_best_fit(pmin(2^LFC_cb,2), Dose)) %>%
    dplyr::ungroup() 
}




# Create a cluster
cl <- makeCluster(detectCores() - 1)
# Fit the curves
DRC <- parLapply(cl, LFC.parallel, f)
# Stop the cluster
stopCluster(cl)

DRC <- dplyr::bind_rows(DRC)
rm(LFC.parallel, f)


DRC <- ViabilityCollapsedConditions %>%
  dplyr::group_by(SampleID) %>% 
  dplyr::summarise(md = min(Dose, na.rm = T),
                   MD = max(Dose, na.rm = T)) %>%
  dplyr::ungroup() %>% 
  dplyr::inner_join(DRC) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(AUC = ifelse(successful_fit, compute_auc(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD), NA),
                log2.IC50 = ifelse(successful_fit, compute_log_ic50(Lower_Limit, Upper_Limit, Inflection, Slope, md, MD), NA)) 



# Writing Files ----

DRC %>% 
  dplyr::rename(ModelID = depmap_id,
                  EC50 = Inflection,
                  LowerAsymptote = Lower_Limit,
                  UpperAsymptote = Upper_Limit) %>% 
  dplyr::distinct(ModelID, SampleID, CompoundPlate, EC50, LowerAsymptote, UpperAsymptote, Slope) %>% 
  write_csv(paste0(path, "secondary-screen-response-curves.csv"))

DRC %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "AUC") %>% 
  write.csv(paste0(path, "secondary-screen-AUC-matrix.csv"))
  
DRC %>% 
  reshape2::acast(depmap_id ~ SampleID, value.var = "log2.IC50") %>% 
  write.csv(paste0(path, "secondary-screen-log2-IC50-matrix.csv"))


ViabilityMatrix %>% 
  write.csv(paste0(path, "secondary-screen-viability-matrix.csv"))

ViabilityCollapsedMatrix %>% 
  write.csv(paste0(path, "secondary-screen-viability-collapsed-matrix.csv"))

ViabilityConditions %>% 
  write_csv(paste0(path, "secondary-screen-viability-conditions.csv"))

ViabilityCollapsedConditions %>% 
  write_csv(paste0(path, "secondary-screen-viability-collapsed-conditions.csv"))
