
#-------------------Header------------------------------------------------
# Author: NAME
# Date: 2/11/21
# Purpose: Generate example plots of model fit
#          
# source("FILEPATH/model_fit_plot.R", echo=T)
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

packages <- c("data.table","magrittr","ggplot2", "matrixStats", "cowplot", "ggpubr", "MASS")

library(mrbrt001, lib.loc = "FILEPATH")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Directories -------------------------------------------------------------

model_name <- "VERSION_NAME"

n_draws <- 1000

work_dir <- file.path(j_root, "FILEPATH", paste0("cumulative_", model_name))
output_dir <- file.path(j_root, "FILEPATH", model_name)
dir.create(output_dir)

source(file.path(central_lib, "FILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)

color_scheme_reg <- c("Central Europe, Eastern Europe, and Central Asia" = "#0A2F51",
                      "High-income" = "#57c1e4",
                      "Latin America and Caribbean" = "#2e9598",
                      "North Africa and Middle East" = "#f7db69",
                      "Southeast Asia, East Asia, and Oceania" = "#f26a44",
                      "South Asia" = "#ec1b4b",
                      "Sub-Saharan Africa" = "#a8216b",
                      "Global" = "grey")

antigens <- c("dtp", "mcv")

# All countries -----------------------------------------------------------

min_c <- 0
max_c <- 1.1

input_data <- fread(paste0(work_dir, "/stage1/input_data.csv"))

pdf(file.path(output_dir, "model_fit_plots_all.pdf"), width = 12, height = 4)

out_pred_loc_range <- data.table()

for(l_id in hierarchy[level == 3][order(super_region_name), unique(location_id)]){
  
  for(ant in antigens){
    
    l_name <- hierarchy[location_id == l_id, location_name]
    sr_name <- hierarchy[location_id == l_id, super_region_name]
    sr_id <- hierarchy[location_id == l_id, super_region_id]
    
    print(paste(l_name, sr_name, ant))
    
    gglist <- list()
    
    # Step 1 input data and prediction (with uncertainty) --------------------
    
    input_data1 <- fread(paste0(work_dir, "/stage1/input_data.csv"))[super_region_id == sr_id & antigen == ant]
    
    input_data1 <- input_data1[, .(location_id, date, month, source_id, source, ratio, log_ratio_se, weight, lower, upper, log_ratio, cum_mob,
                                   super_region_id, super_region_name, model = ant)]
    
    # generate model dataset
    pred1_dt <- data.table(location_id = l_id, cum_mob = seq(0,0.7, 0.01))
    
    dat_pred <- MRData()
    
    dat_pred$load_df(
      data = pred1_dt,
      col_covs=list("cum_mob")
    )
    
    if(l_id %in% input_data1$location_id){
      # predict location specific with uncertainty
      model1 <- py_load_object(paste0(work_dir, "/stage1/pickles/loc_antigen__", l_id,"_", ant, ".pkl"), pickle = "dill")
      
      loc_sample1 <- mrbrt001::core$other_sampling$sample_simple_lme_beta(sample_size = as.integer(n_draws), model = model1)
      
      loc_draws1 <- model1$create_draws(dat_pred, 
                                        beta_samples = loc_sample1, 
                                        gamma_samples = matrix(nrow = n_draws, ncol = 6, data = 0),
                                        random_study = TRUE,                                         
                                        sort_by_data_id = TRUE) %>% exp()
      

      pred1_dt[, location := model1$predict(dat_pred, predict_for_study = F) %>% exp()]
      pred1_dt[, loc_lower := model1$predict(dat_pred, predict_for_study = F) %>% exp()]
      pred1_dt[, loc_upper := model1$predict(dat_pred, predict_for_study = F) %>% exp()]
      
      
      out_pred_loc_range <- rbind(out_pred_loc_range, 
                                  data.table(location_id = l_id,
                                             pred = pred1_dt[cum_mob == 0.4, location],
                                             model = ant)) 
      
    }else{
      
      pred1_dt[, c("loc_lower", "loc_upper", "location") := as.numeric(NA)]
      pred1_dt[cum_mob == 0, c("loc_lower", "loc_upper", "location") := 0]
      
    }
    
    
    # predict region specific
    model1r <- py_load_object(paste0(work_dir, "/stage1/pickles/sr_antigen__", sr_id, "_", ant, ".pkl"), pickle = "dill")
    
    pred1_dt[, region :=  model1r$predict(dat_pred, predict_for_study = F) %>% exp()]
    
    # predict global
    model1g <- py_load_object(paste0(work_dir, "/stage1/pickles/stage1__stage1.pkl"), pickle = "dill")
    
    pred1_dt[, global :=  model1g$predict(dat_pred, predict_for_study = F) %>% exp()]
    
    pred1_dt[, sr_name := sr_name]
    input_data1[, sr_name := sr_name]
    
    final_pred_sum <- fread(paste0(work_dir, "/vaccine_", ant, "/summary_results.csv"))[location_id == l_id] %>% unique()
    input_data3 <- fread(paste0(work_dir, "/input_data_phase_2.csv"))[location_id == l_id & antigen == ant]
    
    final_pred_sum <- rbind(data.table(month = c(1,2), location_id = l_id, location_name = l_name, super_region_name = sr_name,
                                       super_region_id = sr_id, month_ratio_1_mean = 1, month_ratio_2_mean = 1, month_ratio_2_lower = 1,
                                       month_ratio_2_upper = 1),
                            final_pred_sum, fill = T)
    
    max_a <- max(c(pred1_dt$location, input_data1[location_id == l_id, ratio], pred1_dt$global, pred1_dt$region, 
                   final_pred_sum$month_ratio_1_mean, final_pred_sum$month_ratio_2_upper, input_data3$upper_inst))
    min_a <- min(c(pred1_dt$location, input_data1[location_id == l_id, ratio], pred1_dt$global, pred1_dt$region, 
                   final_pred_sum$month_ratio_1_mean, final_pred_sum$month_ratio_2_lower, input_data3$lower_inst))
    
    
    pred1_dt[loc_upper > max_a, loc_upper := Inf]
    
    input_data1[ratio > max_a, ratio := Inf]
    
    pred1_dt[loc_lower < min_a, loc_lower := -Inf]
    
    input_data1[ratio < min_a, ratio := -Inf]
    
    
    gg1 <- ggplot(pred1_dt, aes(x = cum_mob))+
      geom_ribbon(aes(ymin = loc_lower, ymax = loc_upper, fill = sr_name), alpha = 0.2)+
      geom_hline(aes(yintercept = 1), linetype = "dotted")+
      geom_line(data = pred1_dt, aes(y = global, color = "Global"))+
      geom_point(data = input_data1, aes(y = ratio, size = weight, color = sr_name), alpha = 0.2)+
      geom_line(data = pred1_dt, aes(y = region, color = sr_name), alpha = 0.3)+
      geom_line(aes(color = sr_name, y = location), alpha = 0.9)+  # y = loc_mean
      geom_point(data = input_data1[location_id == l_id], aes(y = ratio, size = weight, fill = sr_name), 
                 shape = 21, color = "black", alpha = 0.9)+
      scale_color_manual(values = color_scheme_reg)+
      scale_fill_manual(values = color_scheme_reg)+
      scale_size(range = c(2,5))+
      scale_y_continuous(limits = c(min_a, max_a))+
      scale_x_continuous(expand = c(0,0))+
      labs(title = "A) Step 1",
           x = "Cumulative Mobility Disruption",
           y = "Cumulative Vaccine Disruption Ratio")+
      theme_bw()+
      theme(legend.position = "none")
    
    gglist[[1]] <- gg1
    
    
    # Step 2 input data and prediction ---------------------------------------
    
    input_data2 <- fread(paste0(work_dir, "/resid/input_data.csv"))[super_region_id == sr_id & antigen == ant]
    
    input_data2 <- input_data2[, .(location_id, date, month, source_id, source, logit_resid, logit_ratio_se,
                                   super_region_id, super_region_name, model = ant)]
    
    # generate model dataset
    pred2_dt <- data.table(location_id = l_id, month = 1:12)
    
    dat_pred <- MRData()
    
    dat_pred$load_df(
      data = pred2_dt,
      col_covs=list("month")
    )
    
    
    if(l_id %in% input_data2$location_id){
      # predict location specific with uncertainty
      model2 <- py_load_object(paste0(work_dir, "/resid/pickles/loc_antigen__", l_id, "_", ant, ".pkl"), pickle = "dill")
      
      model2_parent <- py_load_object(file.path(work_dir, "resid/pickles", 
                                                paste0("sr_antigen__", sr_id, "_", ant, ".pkl")), pickle = "dill")
      
      model2_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model = model2)
      vcov2 <- solve(mrbrt001::core$other_sampling$extract_simple_lme_hessian(model2_specs))
      
      model2_parent_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model = model2_parent)
      vcov2_parent <- solve(mrbrt001::core$other_sampling$extract_simple_lme_hessian(model2_parent_specs))
      
      scalar2 <- max(diag(vcov2)/diag(vcov2_parent))
      
      loc_sample_2 <- mvrnorm(n = n_draws, mu = matrix(model2$beta_soln), Sigma = vcov2/scalar2)
      
      loc_draws2 <- model2$create_draws(dat_pred, 
                                        beta_samples = loc_sample_2, 
                                        gamma_samples = matrix(nrow = n_draws, ncol = 1, data = 0),
                                        random_study = TRUE,                                         
                                        sort_by_data_id = TRUE)
      
      pred2_dt[, location := rowMeans(loc_draws2)]
      pred2_dt[, loc_lower := rowQuantiles(loc_draws2, prob = 0.025)]
      pred2_dt[, loc_upper := rowQuantiles(loc_draws2, prob = 0.975)]
      
    }else{
      pred2_dt[, c("loc_lower", "loc_upper", "location") := as.numeric(NA)]
      pred2_dt[month == 0, c("loc_lower", "loc_upper", "location") := 0]
    }
    
    
    
    
    # predict region specific
    model2r <- py_load_object(paste0(work_dir, "/resid/pickles/sr_antigen__", sr_id, "_", ant, ".pkl"), pickle = "dill")
    
    pred2_dt[, region :=  model2r$predict(dat_pred, predict_for_study = F)]
    
    # predict global
    model2g <- py_load_object(paste0(work_dir,  "/resid/pickles/stage1__stage1.pkl"), pickle = "dill")
    
    pred2_dt[, global :=  model2g$predict(dat_pred, predict_for_study = F)]
    
    pred2_dt[, sr_name := sr_name]
    input_data2[, sr_name := sr_name]
    
    gg2 <- ggplot(pred2_dt, aes(x = month, y = location))+
      geom_ribbon(aes(ymin = loc_lower, ymax = loc_upper, fill = sr_name), alpha = 0.2)+
      geom_hline(aes(yintercept = 0), linetype = "dotted")+
      geom_line(aes(y = global, color = "Global"))+
      geom_point(data = input_data2, aes(y = logit_resid, size = logit_ratio_se^(-2), color = sr_name), alpha = 0.2)+
      geom_line(aes(y = region, color = sr_name), alpha = 0.3)+
      geom_line(aes(color = sr_name), alpha = 0.9)+
      geom_point(data = input_data2[location_id == l_id], aes(y = logit_resid, fill = sr_name, size = logit_ratio_se^(-2)), 
                 shape = 21, color = "black", alpha = 0.9)+
      scale_color_manual(values = color_scheme_reg)+
      scale_fill_manual(values = color_scheme_reg)+
      scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                         expand = c(0,0))+
      scale_size(range = c(2,5))+
      labs(title = "B) Step 2",
           x = "Month",
           y = "Residual (logit space)")+
      theme_bw()+
      theme(legend.position = "none")
    
    gglist[[2]] <- gg2
    
    
    # Mobility ----------------------------------------------------------------
    
    mob <- fread(paste0(work_dir, "/pred_frame.csv"))[location_id == l_id]
    
    mob[, mob_avg_month := 1-mob_avg_month]
    
    mob[mob_avg_month > max_c, mob_avg_month := Inf]
    
    mob[mob_avg_month < min_c, mob_avg_month := -Inf]
    
    mob[, sr_name := sr_name]
    
    gg3 <- ggplot(data = mob, aes(x = month))+
      geom_hline(aes(yintercept = 1), linetype = "dotted")+
      geom_line(aes(y = mob_avg_month, color = sr_name))+
      scale_color_manual(values = color_scheme_reg)+
      scale_y_continuous(limits = c(min_c, max_c))+
      labs(title = "C) Monthly Mobility",
           y = "Monthly Vaccine Disruption Ratio",
           x = "Month")+
      scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                         expand = c(0,0))+
      theme_bw()+
      theme(legend.position = "none")
    
    gglist[[3]] <- gg3
    
    # Model prediction monthly space ------------------------------------------
    
    
    final_pred_sum[, sr_name := sr_name]
    input_data3[, sr_name := sr_name]
    
    max_d <- max_a
    min_d <- min_a
    
    final_pred_sum[month_ratio_2_upper > max_d, month_ratio_2_upper := Inf]
    
    input_data3[inst_ratio > max_d, ratio := Inf]
    input_data3[lower_inst > max_d, ratio := Inf]
    input_data3[upper_inst > max_d, ratio := Inf]
    
    final_pred_sum[month_ratio_2_lower < min_d, month_ratio_2_lower := -Inf]
    
    input_data3[inst_ratio < min_d, ratio := -Inf]
    input_data3[lower_inst < min_d, ratio := -Inf]
    input_data3[upper_inst < min_d, ratio := -Inf]
    
    gg4 <- ggplot(data = final_pred_sum, aes(x = month))+
      geom_hline(aes(yintercept = 1), linetype = "dotted")+
      geom_line(aes(y = month_ratio_1_mean, linetype = "2", color = sr_name))+
      geom_line(aes(y = month_ratio_2_mean, linetype = "1", color = sr_name))+
      geom_ribbon(aes(ymin = month_ratio_2_lower, ymax = month_ratio_2_upper, fill = sr_name), alpha = 0.2)+
      geom_point(data = input_data3, aes(y = inst_ratio, color = sr_name), alpha = 0.9)+
      geom_errorbar(data = input_data3, aes(ymin = lower_inst, ymax = upper_inst, color = sr_name), alpha = 0.9)+
      scale_color_manual(values = color_scheme_reg)+
      scale_fill_manual(values = color_scheme_reg)+
      scale_y_continuous(limits = c(min_d, max_d))+
      labs(title = "D) Monthly Results",
           y = "Monthly Vaccine Disruption Ratio",
           x = "Month")+
      scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                         expand = c(0,0))+
      theme_bw()+
      theme(legend.position = "none")
    
    gglist[[4]] <- gg4
    
    p_all <- plot_grid(plotlist = gglist, ncol = 4)
    
    main_title <- ggdraw() + draw_label(paste0(l_name, ", ", toupper(ant), ", ", sr_name), fontface='bold')
    
    title <- plot_grid(main_title, ncol = 1)
    
    p_all <- plot_grid(title, p_all, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
    
    print(p_all)
    
  }
  
}


dev.off()

out_pred_loc_range[, .(min(pred), max(pred), sd(pred)), by = "model"]



# Paper locs --------------------------------------------------------------


locs <- c("dtp" = 163, "mcv"= 102, "dtp" = 214, "mcv" = 164, "dtp" = 71)

min_as <- c(0.35, 0.4, 0.6, 0.2, 0.7)
max_as <- c(1.05, 1.05, 1.05, 1.5, 1.05)


min_b <- -5
max_b <- 5

min_c <- 0
max_c <- 1.1


gglist <- list()

for(i in 1:length(locs)){
  
  min_a <- min_as[i]
  max_a <- max_as[i]
  
  min_d <- min_a
  max_d <- max_a
  
  l_id <- locs[i]
  vacc <- names(locs)[i]
  ant <- vacc
  
  if(vacc == "dtp"){
    antigen <- "DTP3"
  }else if(vacc == "mcv"){
    antigen <- "MCV1"
  }
  
  l_name <- hierarchy[location_id == l_id, location_name]
  sr_name <- hierarchy[location_id == l_id, super_region_name]
  sr_id <- hierarchy[location_id == l_id, super_region_id]
  
  print(paste(l_name, sr_name, vacc))
  
  # Step 1 input data and prediction  --------------------
  
  input_data1 <- fread(paste0(work_dir, "/stage1/input_data.csv"))[super_region_id == sr_id & antigen == vacc]
  
  input_data1 <- input_data1[, .(location_id, date, month, source_id, source, ratio, log_ratio_se, weight, lower, upper, log_ratio, cum_mob,
                                 super_region_id, super_region_name, model = vacc)]
  
  # generate model dataset
  pred1_dt <- data.table(location_id = l_id, cum_mob = seq(0,0.7, 0.01))
  
  
  dat_pred <- MRData()
  
  dat_pred$load_df(
    data = pred1_dt,
    col_covs=list("cum_mob")
  )
  
  if(l_id %in% input_data1$location_id){
    # predict location specific with uncertainty
    model1 <- py_load_object(paste0(work_dir, "/stage1/pickles/loc_antigen__", l_id,"_", ant, ".pkl"), pickle = "dill")
    
    loc_sample1 <- mrbrt001::core$other_sampling$sample_simple_lme_beta(sample_size = as.integer(n_draws), model = model1)
    
    loc_draws1 <- model1$create_draws(dat_pred, 
                                      beta_samples = loc_sample1, 
                                      gamma_samples = matrix(nrow = n_draws, ncol = 6, data = 0),
                                      random_study = TRUE,                                         
                                      sort_by_data_id = TRUE) %>% exp()
    
    pred1_dt[, location := model1$predict(dat_pred, predict_for_study = F) %>% exp()]
    pred1_dt[, loc_lower := model1$predict(dat_pred, predict_for_study = F) %>% exp()]
    pred1_dt[, loc_upper := model1$predict(dat_pred, predict_for_study = F) %>% exp()]
    
    
  }else{
    
    pred1_dt[, c("loc_lower", "loc_upper", "location") := as.numeric(NA)]
    pred1_dt[cum_mob == 0, c("loc_lower", "loc_upper", "location") := 0]
    
  }
  
  
  # predict region specific
  model1r <- py_load_object(paste0(work_dir, "/stage1/pickles/sr_antigen__", sr_id, "_", vacc, ".pkl"), pickle = "dill")
  
  pred1_dt[, region :=  model1r$predict(dat_pred, predict_for_study = F) %>% exp()]
  
  # predict global
  model1g <- py_load_object(paste0(work_dir, "/stage1/pickles/stage1__stage1.pkl"), pickle = "dill")
  
  pred1_dt[, global :=  model1g$predict(dat_pred, predict_for_study = F) %>% exp()]
  
  pred1_dt[, sr_name := sr_name]
  input_data1[, sr_name := sr_name]
  
  pred1_dt[loc_upper > max_a, loc_upper := Inf]
  
  input_data1[ratio > max_a, ratio := Inf]
  
  pred1_dt[loc_lower < min_a, loc_lower := -Inf]
  
  input_data1[ratio < min_a, ratio := -Inf]
  
  
  gg1 <- ggplot(pred1_dt, aes(x = cum_mob))+
    geom_hline(aes(yintercept = 1), linetype = "dotted")+
    geom_line(data = pred1_dt, aes(y = global, color = "Global"))+
    geom_point(data = input_data1, aes(y = ratio, size = weight, color = sr_name), alpha = 0.2)+
    geom_line(data = pred1_dt, aes(y = region, color = sr_name), alpha = 0.3)+
    geom_line(aes(color = sr_name, y = location), alpha = 0.9)+  # y = loc_mean
    geom_point(data = input_data1[location_id == l_id], aes(y = ratio, size = weight, fill = sr_name),
               shape = 21, color = "black", alpha = 0.9)+
    scale_color_manual(values = color_scheme_reg)+
    scale_fill_manual(values = color_scheme_reg)+
    scale_size(range = c(2,5))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(limits = c(min_a, max_a), 
                       breaks = seq(0,2, 0.2))+
    labs(subtitle = paste0(l_name, ", ", antigen),
         x = "",
         y = "")+
    theme_bw()+
    theme(legend.position = "none")
  
  gglist[[(4*i-3)]] <- gg1
  
  
  # Step 2 input data and prediction ---------------------------------------
  
  input_data2 <- fread(paste0(work_dir, "/resid/input_data.csv"))[super_region_id == sr_id & antigen == vacc]
  
  input_data2 <- input_data2[, .(location_id, date, month, source_id, source, logit_resid, logit_ratio_se,
                                 super_region_id, super_region_name, model = vacc)]
  
  # generate model dataset
  pred2_dt <- data.table(location_id = l_id, month = 1:12)
  
  dat_pred <- MRData()
  
  dat_pred$load_df(
    data = pred2_dt,
    col_covs=list("month")
  )
  
  
  if(l_id %in% input_data2$location_id){
    # predict location specific with uncertainty
    model2 <- py_load_object(paste0(work_dir, "/resid/pickles/loc_antigen__", l_id, "_", ant, ".pkl"), pickle = "dill")
    
    model2_parent <- py_load_object(file.path(work_dir, "resid/pickles", 
                                              paste0("sr_antigen__", sr_id, "_", ant, ".pkl")), pickle = "dill")
    
    model2_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model = model2)
    vcov2 <- solve(mrbrt001::core$other_sampling$extract_simple_lme_hessian(model2_specs))
    
    model2_parent_specs <- mrbrt001::core$other_sampling$extract_simple_lme_specs(model = model2_parent)
    vcov2_parent <- solve(mrbrt001::core$other_sampling$extract_simple_lme_hessian(model2_parent_specs))
    
    scalar2 <- max(diag(vcov2)/diag(vcov2_parent))
    
    loc_sample_2 <- mvrnorm(n = n_draws, mu = matrix(model2$beta_soln), Sigma = vcov2/scalar2)
    
    loc_draws2 <- model2$create_draws(dat_pred, 
                                      beta_samples = loc_sample_2, 
                                      gamma_samples = matrix(nrow = n_draws, ncol = 1, data = 0),
                                      random_study = TRUE,                                         
                                      sort_by_data_id = TRUE)
    
    pred2_dt[, location := rowMeans(loc_draws2)]
    pred2_dt[, loc_lower := rowQuantiles(loc_draws2, prob = 0.025)]
    pred2_dt[, loc_upper := rowQuantiles(loc_draws2, prob = 0.975)]
    
  }else{
    pred2_dt[, c("loc_lower", "loc_upper", "location") := as.numeric(NA)]
    pred2_dt[month == 0, c("loc_lower", "loc_upper", "location") := 0]
  }
  
  # predict region specific
  model2r <- py_load_object(paste0(work_dir, "/resid/pickles/sr_antigen__", sr_id, "_", vacc, ".pkl"), pickle = "dill")
  
  pred2_dt[, region :=  model2r$predict(dat_pred, predict_for_study = F)]
  
  # predict global
  model2g <- py_load_object(paste0(work_dir, "/resid/pickles/stage1__stage1.pkl"), pickle = "dill")
  
  pred2_dt[, global :=  model2g$predict(dat_pred, predict_for_study = F)]
  
  pred2_dt[, sr_name := sr_name]
  input_data2[, sr_name := sr_name]
  
  input_data2[logit_resid > max_b, logit_resid := Inf]
  
  input_data2[logit_resid < min_b, logit_resid := -Inf]
  
  gg2 <- ggplot(pred2_dt, aes(x = month, y = location))+
    geom_hline(aes(yintercept = 0), linetype = "dotted")+
    geom_line(aes(y = global, color = "Global"))+
    geom_point(data = input_data2, aes(y = logit_resid, size = logit_ratio_se^(-2), color = sr_name), alpha = 0.2)+
    geom_line(aes(y = region, color = sr_name), alpha = 0.3)+
    geom_line(aes(color = sr_name), alpha = 0.9)+
    geom_point(data = input_data2[location_id == l_id], aes(y = logit_resid, fill = sr_name, size = logit_ratio_se^(-2)),
               shape = 21, color = "black", alpha = 0.9)+
    scale_color_manual(values = color_scheme_reg)+
    scale_fill_manual(values = color_scheme_reg)+
    scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                       expand = c(0,0), minor_breaks = NULL)+
    scale_y_continuous(limits = c(min_b, max_b),
                       breaks = c(-4, -2, 0, 2, 4))+
    scale_size(range = c(2,5))+
    labs(subtitle = " ",
         x = "",
         y = "")+
    theme_bw()+
    theme(legend.position = "none")
  
  gglist[[(4*i-2)]] <- gg2
  
  
  # Mobility ----------------------------------------------------------------
  
  mob <- fread(paste0(work_dir, "/pred_frame.csv"))[location_id == l_id]
  
  mob[, mob_avg_month := 1-mob_avg_month]
  
  mob[mob_avg_month > max_c, mob_avg_month := Inf]
  
  mob[mob_avg_month < min_c, mob_avg_month := -Inf]
  
  mob[, sr_name := sr_name]
  
  gg3 <- ggplot(data = mob, aes(x = month))+
    geom_hline(aes(yintercept = 1), linetype = "dotted")+
    geom_line(aes(y = mob_avg_month, color = sr_name))+
    scale_color_manual(values = color_scheme_reg)+
    scale_y_continuous(limits = c(min_c, max_c),
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1),
                       expand = c(0,0))+
    labs(subtitle = " ",y = "",
         x = "")+
    scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                       expand = c(0,0), minor_breaks = NULL)+
    theme_bw()+
    theme(legend.position = "none")
  
  gglist[[(4*i-1)]] <- gg3
  
  # Model prediction monthly space ------------------------------------------
  
  final_pred_sum <- fread(paste0(work_dir, "/vaccine_", vacc, "/summary_results.csv"))[location_id == l_id] %>% unique()
  input_data3 <- fread(paste0(work_dir, "/input_data_phase_2.csv"))[location_id == l_id & antigen == vacc]
  
  final_pred_sum <- rbind(data.table(month = c(1,2), location_id = l_id, location_name = l_name, super_region_name = sr_name,
                                     super_region_id = sr_id, month_ratio_1_mean = 1, month_ratio_2_mean = 1, month_ratio_2_lower = 1,
                                     month_ratio_2_upper = 1),
                          final_pred_sum, fill = T)
  
  final_pred_sum[, sr_name := sr_name]
  input_data3[, sr_name := sr_name]
  
  final_pred_sum[month_ratio_1_mean > max_d, month_ratio_1_mean := Inf]
  final_pred_sum[month_ratio_2_mean > max_d, month_ratio_2_mean := Inf]
  final_pred_sum[month_ratio_2_lower > max_d, month_ratio_2_lower := Inf]
  final_pred_sum[month_ratio_2_upper > max_d, month_ratio_2_upper := Inf]
  
  input_data3[inst_ratio > max_d, ratio := Inf]
  input_data3[lower_inst > max_d, ratio := Inf]
  input_data3[upper_inst > max_d, ratio := Inf]
  
  final_pred_sum[month_ratio_1_mean < min_d, month_ratio_1_mean := -Inf]
  final_pred_sum[month_ratio_2_mean < min_d, month_ratio_2_mean := -Inf]
  final_pred_sum[month_ratio_2_lower < min_d, month_ratio_2_lower := -Inf]
  final_pred_sum[month_ratio_2_upper < min_d, month_ratio_2_upper := -Inf]
  
  input_data3[inst_ratio < min_d, ratio := -Inf]
  input_data3[lower_inst < min_d, ratio := -Inf]
  input_data3[upper_inst < min_d, ratio := -Inf]
  
  gg4 <- ggplot(data = final_pred_sum, aes(x = month))+
    geom_hline(aes(yintercept = 1), linetype = "dotted")+
    geom_line(aes(y = month_ratio_1_mean, linetype = "2", color = sr_name))+
    geom_line(aes(y = month_ratio_2_mean, linetype = "1", color = sr_name))+
    geom_ribbon(aes(ymin = month_ratio_2_lower, ymax = month_ratio_2_upper, fill = sr_name), alpha = 0.2)+
    geom_point(data = input_data3, aes(y = inst_ratio, color = sr_name), alpha = 0.9)+
    geom_errorbar(data = input_data3, aes(ymin = lower_inst, ymax = upper_inst, color = sr_name), alpha = 0.9)+
    scale_color_manual(values = color_scheme_reg)+
    scale_fill_manual(values = color_scheme_reg)+
    scale_y_continuous(limits = c(min_d, max_d),
                       breaks = seq(0,2, 0.2))+
    labs(subtitle = " ",y = "",
         x = "")+
    scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"),
                       expand = c(0,0), minor_breaks = NULL)+
    theme_bw()+
    theme(legend.position = "none")
  
  gglist[[(4*i)]] <- gg4
  
}


# Build titles and legends ------------------------------------------------

gglist[[9]] <- gglist[[9]] + labs(y = "Cumulative Vaccine Disruption Ratio")
gglist[[10]] <- gglist[[10]] + labs(y = "Residual (logit space)")
gglist[[11]] <- gglist[[11]] + labs(y = "Monthly Mobility Disruption Ratio")
gglist[[12]] <- gglist[[12]] + labs(y = "Monthly Vaccine Disruption Ratio")

gglist[[17]] <- gglist[[17]] + labs(x = "Cumulative Mobility Disruption")
gglist[[18]] <- gglist[[18]] + labs(x = "Month")
gglist[[19]] <- gglist[[19]] + labs(x = "Month")
gglist[[20]] <- gglist[[20]] + labs(x = "Month")

plots <- plot_grid(plotlist = gglist, ncol = 4)

A <- ggdraw() + draw_label("A) Step 1", fontface='bold')
B <- ggdraw() + draw_label("B) Step 2", fontface='bold')
C <- ggdraw() + draw_label("C) Monthly Mobility", fontface='bold')
D <- ggdraw() + draw_label("D) Monthly Results", fontface='bold')

title <- plot_grid(A, B, C, D, ncol = 4)

legend_1 <- 
  ggplot(data.table(a = c("Global, both antigens"), b = c("Country data", "Super-region data")), 
         aes(x = a, y = a))+
  geom_line(aes(color = a))+
  geom_point(aes(shape = b, fill = b, alpha = b), fill = "grey")+
  scale_shape_manual(values = c("Country data" = 21, "Super-region data" = 19))+
  scale_alpha_manual(values = c("Country data" = 0.9, "Super-region data" = 0.2))+
  scale_color_manual(values = c("Global, both antigens" = "grey"))+
  labs(color = "", shape = "", fill = "", alpha = "")+
  theme_bw()+
  theme(legend.spacing.y = unit(0, "cm"))

legend_1 <- get_legend(legend_1)

legend_2_dt <- expand.grid(a = c("High-income", "South Asia", "Sub-Saharan Africa"), 
                           b = c("country-antigen-model", "super-region-antigen-model")) %>% as.data.table()

legend_2_dt[, c:= paste0(a, ", ", b)]

color_scheme_reg_2 <- c("High-income, country-antigen-model" = "#57c1e4",
                        "South Asia, country-antigen-model" = "#ec1b4b",
                        "Sub-Saharan Africa, country-antigen-model" = "#a8216b",
                        "High-income, super-region-antigen-model" = "#57c1e4",
                        "South Asia, super-region-antigen-model" = "#ec1b4b",
                        "Sub-Saharan Africa, super-region-antigen-model" = "#a8216b")

alpha_2 <- c("High-income, super-region-antigen-model" = 0.3,
             "South Asia, super-region-antigen-model" = 0.3,
             "Sub-Saharan Africa, super-region-antigen-model" = 0.3,
             "High-income, country-antigen-model" = 0.9,
             "South Asia, country-antigen-model" = 0.9,
             "Sub-Saharan Africa, country-antigen-model" = 0.9)

legend_2_dt[, c:= factor(c, levels = names(alpha_2))]

legend_2 <- 
  ggplot(legend_2_dt, aes(x = a, y = a))+
  geom_line(aes(color = c, alpha = c))+
  scale_alpha_manual(values = alpha_2)+
  scale_color_manual(values = color_scheme_reg_2)+
  labs(color = "", alpha = "")+
  guides(color=guide_legend(ncol=2), alpha=guide_legend(ncol=2))+
  theme_bw()

legend_2 <- get_legend(legend_2)

legend_3_dt <- data.table(a = c("Step 1", "Model estimates"), b = "Monthly country data with 95% UI")

legend_3 <- 
  ggplot(legend_3_dt, aes(x = a, y = b, ymin = b, ymax = b))+
  geom_line(aes(linetype = a))+
  geom_point(aes(color = b))+
  geom_pointrange(aes(color = b))+
  scale_color_manual(values = "black")+
  scale_linetype(guide = guide_legend(reverse = T))+
  labs(linetype = "", color = "")+
  theme_bw()+
  theme(legend.spacing.y = unit(0, "cm"))

legend_3 <- get_legend(legend_3)

legends <- plot_grid(legend_1, legend_2, legend_3, ncol = 3, rel_widths = c(1,2,1))


# Put it all together -----------------------------------------------------

p_all <- plot_grid(plotlist = list(title, plots, legends), ncol=1, rel_heights=c(0.02, 1, 0.15)) # rel_heights values control title margins


pdf(file.path(output_dir, "model_fit_plots.pdf"), width = 12, height = 12)

print(p_all)

dev.off()


out_pred_loc_range[, .(min(pred), max(pred), sd(pred)), by = "model"]
