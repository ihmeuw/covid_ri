
#-------------------Header------------------------------------------------
# Author: NAME
# Date: 2/3/21
# Purpose: Child script to generate draws from model
#          
# source("FILEPATH/02b_spline_cascade_draws.R", echo=T)
#***************************************************************************

#------------------SET-UP--------------------------------------------------

# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  j_root <- "FILEPATH"
  h_root <- "FILEPATH"
  central_lib <- "FILEPATH"
  } else {
  j_root <- "FILEPATH"
  h_root <- "FILEPATH"
  central_lib <- "FILEPATH"
  }



# load packages, install if missing

lib.loc <- paste0(h_root,"R/",R.Version()$platform,"/",R.Version()$major,".",R.Version()$minor)
dir.create(lib.loc,recursive=T, showWarnings = F)
.libPaths(c(lib.loc,.libPaths()))

packages <- c("data.table","magrittr","matrixStats", "MASS")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

args <- commandArgs(trailingOnly = TRUE)

# set up defaults for running as its own script, not in parallel

if(length(args)==0){
  args <- list(100,522, FALSE,
               "FILEPATH",
               "dtp", 0.001, 10)
}


library(mrbrt001, lib.loc = "FILEPATH")

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

# Directories -------------------------------------------------------------

# arguments
n_draws <- args[[1]] %>% as.integer()
loc_id <- args[[2]] %>% as.integer()
super_region <- args[[3]] %>% as.logical()
work_dir <- args[[4]]
ant <- args[[5]]
offset <- args[[6]] %>% as.numeric
region_weight <- args[[7]] %>% as.numeric

set.seed(145)

ratio <- paste0("vaccine_", ant)

draw_cols <- paste0("draw_", 0:(n_draws-1))

final_pred_dt <- fread(file.path(work_dir, "pred_frame.csv"))
dt <- fread(file.path(work_dir, "input_data.csv"))

loc_sample_1 <- matrix(nrow = n_draws, ncol = 6)
loc_sample_2 <- matrix(nrow = n_draws, ncol = 4)

if(super_region){
  
  tmp_pred <- final_pred_dt[!(location_id %in% dt$location_id) & super_region_id == loc_id & month %in% c(3:12)]
  
  # sample from country-specific regional and global results
  sample_locs <- sample(c(unique(dt[, location_id]), rep(unique(dt[super_region_id == loc_id, location_id]), region_weight-1)), size = n_draws, replace = T)
  
  i <- 1
  
  for(loc in unique(sample_locs)){
    
    n_samples <- sum(sample_locs == loc)
    
    model1 <- py_load_object(file.path(work_dir, "stage1/pickles", paste0("loc_antigen__", loc, "_", ant,".pkl")), pickle = "dill")
    
    loc_sample_1[i:(i+n_samples-1),] <- mrbrt001::core$other_sampling$sample_simple_lme_beta(sample_size = as.integer(n_samples), model = model1)
       
    i <- i+n_samples
    
  }
  
  if(loc_id == 31){
    model2 <- py_load_object(file.path(work_dir, "resid/pickles", paste0("stage1__stage1.pkl")), pickle = "dill")
  }else{
    model2 <- py_load_object(file.path(work_dir, "resid/pickles", paste0("sr_antigen__", loc_id, "_", ant, ".pkl")), pickle = "dill")
  }
  
  loc_sample_2 <- mrbrt001::core$other_sampling$sample_simple_lme_beta(sample_size = as.integer(n_draws), model = model2)
 
  
  
}else{
  
  tmp_pred <- final_pred_dt[location_id == loc_id & month %in% c(3:12)]
  
  # sample directly from location-specific results
  model1 <- py_load_object(file.path(work_dir, "stage1/pickles", paste0("loc_antigen__", loc_id, "_", ant,".pkl")), pickle = "dill")
  
  loc_sample_1 <- mrbrt001::core$other_sampling$sample_simple_lme_beta(sample_size = as.integer(n_draws), model = model1)
  
  if(loc_id == 92){
    model2 <- py_load_object(file.path(work_dir, "resid/pickles", paste0("sr_antigen__64_", ant, ".pkl")), pickle = "dill")
  }else{
    model2 <- py_load_object(file.path(work_dir, "resid/pickles", paste0("loc_antigen__", loc_id, "_", ant,".pkl")), pickle = "dill")
  }
  
  model2_parent <- py_load_object(file.path(work_dir, "resid/pickles", 
                                            paste0("sr_antigen__", unique(tmp_pred$super_region_id), "_", ant, ".pkl")), pickle = "dill")
  
  model2_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model = model2)
  vcov2 <- solve(mrbrt001::core$other_sampling$extract_simple_lme_hessian(model2_specs))
  
  model2_parent_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model = model2_parent)
  vcov2_parent <- solve(mrbrt001::core$other_sampling$extract_simple_lme_hessian(model2_parent_specs))
  
  scalar2 <- max(diag(vcov2)/diag(vcov2_parent))
  
  loc_sample_2 <- mvrnorm(n = n_draws, mu = matrix(model2$beta_soln), Sigma = vcov2/scalar2)
  
}



  dat_pred <- MRData()
  
  dat_pred$load_df(
    data = tmp_pred,
    col_covs=list("cum_mob", "month")
  )
  
  loc_draws <- model1$create_draws(dat_pred, 
                                   beta_samples = loc_sample_1, 
                                   gamma_samples = matrix(nrow = n_draws, ncol = 6, data = 0),
                                   random_study = TRUE,                                         
                                   sort_by_data_id = TRUE) %>% exp()
  
  loc_draws_2 <- model2$create_draws(dat_pred, 
                                     beta_samples = loc_sample_2, 
                                     gamma_samples = matrix(nrow = n_draws, ncol = 1, data = 0),
                                     random_study = TRUE,                                         
                                     sort_by_data_id = TRUE)
  
  tmp_pred[, c(paste0("cum_", draw_cols)) := as.data.table(loc_draws)]
  
  tmp_pred[, c(paste0("resid_", draw_cols)) := as.data.table(loc_draws_2)]
  
  write.csv(tmp_pred, file.path(work_dir, ratio, "draws", paste0(loc_id, "_draws.csv")), row.names = F)
  
  
  


