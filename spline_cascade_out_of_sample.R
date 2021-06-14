#-------------------Header------------------------------------------------
# Author: NAME
# Date: 5/6/2021
# Purpose: Run out of sample validation on spline cascade
#          
# source("FILEPATH/spline_cascade_out_of_sample.R", echo=T)
#***************************************************************************

#------------------SET-UP--------------------------------------------------

# clear memory
rm(list=ls())

start_time <- Sys.time()

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

project <- "-P proj_sdg " 
sge.output.dir <- " -o FILEPATH -e FILEPATH "
rshell <- "FILEPATH/execRscript.sh"


# load packages, install if missing

lib.loc <- paste0(h_root,"R/",R.Version()$platform,"/",R.Version()$major,".",R.Version()$minor)
dir.create(lib.loc,recursive=T, showWarnings = F)
.libPaths(c(lib.loc,.libPaths()))

packages <- c("data.table","magrittr","ggplot2", "msm", "dplyr", "gridExtra", "ggpubr", "wCorr", 
              "dplyr", "matrixStats", "scales", "timeDate", "cowplot")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

library(mrbrt001, lib.loc = "FILEPATH")

# Directories -------------------------------------------------------------

mobility_version <- "VERSION_NUBMVER"  # Change this to reflect the correct version of mobility from Steve and the covid team 

data_dir <-  "VERSION_NUMBER" #Sys.Date() #

model_name <- "MODEL_NAME"

thetas_1 <- c(4, 1, 100)
thetas_2 <- c(20, 1, 100)


# CHANGE THESE IF YOU NEED UNCERTAINTY
uncertainty <- TRUE 
n_draws <- 1000  # Number of draws for calculating uncertainty (only used if uncertainty is true)
draw_cols <- paste0("draw_", 0:(n_draws - 1))
region_weight <- 10 # for countries with no data, countries from the same region are X times as likely to be sampled as all other countries

set.seed(145)

offset <- 1e-3

last_days <- as.Date(timeLastDayInMonth(paste0("2020-", 1:12, "-01")))

end_week_last_days <-  as.Date(c("2020-02-01", "2020-02-29", "2020-03-28", "2020-05-02", "2020-05-30", "2020-06-27", 
                                 "2020-08-01", "2020-08-29", "2020-10-03", "2020-10-31", "2020-11-28", "2021-01-02"))

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

work_root <- paste0(j_root, "FILEPATH", "cumulative_", Sys.Date(), "_", model_name)  # Outputs will be stored here
dir.create(work_root)

vaccine_data_path <- file.path(j_root, "FILEPATH", data_dir, "data_mrbrt_source.csv")

source("FILEPATH/utils.R")

source(file.path(central_lib, "cFILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)

countries <-hierarchy[level == 3, location_id]


# #get population for US, india, and UK subnats
pop <- fread("FILEPATH/all_populations.csv")
pop <- pop[sex_id == 3 & age_group_id == 22, 
           .(location_id, population)]

# Prepare cumulative mobility for prediciton ---------------------------------

mob <- fread(file.path("FILEPATH", mobility_version, "mobility/reference_scenario.csv"))
mob <- mob[, .(location_id, date = as.Date(date), mob_avg = mobility_reference * -0.01)]

# Set Jan and Feb to be zero except China and Taiwan
mob[date < as.Date("2020-03-01") & !(location_id %in% c(6,8)), mob_avg := 0]

# cap mobility at zero disruption
mob[mob_avg<0, mob_avg := 0]

# create cumulative mobility
mob <- mob[order(location_id, date)]
mob[, cum_mob:= 0]
mob <- mob[date>="2020-03-01", cum_mob := cummean(mob_avg), by = "location_id"]
mob[location_id %in% c(6,8), cum_mob := cummean(mob_avg), by = "location_id"]

mob_nat <- mob[location_id %in% countries & date < as.Date("2021-01-01")]

mob_nat[, month:= month(date)]
mob_nat[, mob_avg_month := mean(mob_avg), by = c("month", "location_id")]


# Read in data and split countries into groups ----------------------------

dt_all <- fread(vaccine_data_path)

dt_all[, date := as.Date(end_date)]

# mobility
dt_all <- merge(dt_all, mob, by = c("location_id", "date"), all.x = T)

# Change all subnat locations to parent
for(loc_id in unique(dt_all$location_id)){
  if(!loc_id %in% hierarchy[level ==3, location_id]){
    
    iso3 <- strsplit(dt_all[location_id == loc_id, unique(ihme_loc_id)], "_")[[1]][1]
    
    dt_all[location_id == loc_id, location_id:= hierarchy[ihme_loc_id == iso3, location_id]]
  }
}

dt_all <- merge(dt_all, hierarchy[,.(location_id, super_region_id, super_region_name)], by = "location_id", all.x = T)

write.csv(dt_all, file.path(work_root, "input_data.csv"), row.names = F)

# create random list of holdouts

loc_list <- unique(dt_all[, .(location_id, super_region_id)])

loc_list[, rand:= runif(nrow(loc_list))]

loc_list <- loc_list[order(super_region_id, rand)]

loc_list[, holdout := rep(1:5, ceiling(nrow(loc_list)/5))[1:nrow(loc_list)]]

write.csv(loc_list, file.path(work_root, "holdout_locs.csv"), row.names = F)

# generat count of months of data working backwards
dt_all[order(-month), last_month_n := 1:.N, by = c("ihme_loc_id", "dtp", "mcv", "source_id")]

# Run MRBRT ---------------------------------------------------------------

for(n in 1:10){ # number of holdout months
  for(i in 1:5){ # randomized 20% holdouts
    
    work_dir <- file.path(work_root, paste0("months_", n), paste0("holdout_", i))
    dir.create(work_dir, recursive = T)
    
    pdf(paste0(work_dir, "/spline_fit_plots.pdf"), width = 11, height = 8.5)
    
    
    dt <- dt_all[!(location_id %in% loc_list[holdout == i, location_id] & last_month_n <= n)]
    
    dtp <- dt[dtp == 1]
    dtp[, antigen := "dtp"]
    
    mcv <- dt[mcv == 1]
    mcv[, antigen := "mcv"]
    
    dt <- rbind(dtp, mcv)
    
    setnames(dt, "ratio_se_log", "log_ratio_se")
    
    dt[, loc_antigen := paste0(location_id, "_", antigen)]
    dt[, sr_antigen := paste0(super_region_id, "_", antigen)]
    
    dat_1 <- MRData()
    dat_1$load_df(
      data = dt,  col_obs = "log_ratio", col_obs_se = "log_ratio_se",
      col_covs = list("cum_mob"), 
      col_study_id = "loc_antigen"
    )
    
    max_mob <- max(dt$cum_mob)
    
    best_knots1 <- c(0, quantile(dt$cum_mob, c(1/5, 2/5, 3/5, 4/5))/max_mob, 1/max_mob)
    
    
    model <- MRBRT(data = dat_1,
                   cov_models = list(LinearCovModel("cum_mob",
                                                    use_spline = TRUE,
                                                    use_re = FALSE,
                                                    spline_degree = 3L,
                                                    prior_spline_monotonicity = "decreasing",
                                                    spline_knots_type = "domain",
                                                    spline_knots = array(best_knots1*max_mob),
                                                    spline_r_linear = TRUE,
                                                    spline_l_linear = FALSE,
                                                    prior_spline_maxder_gaussian = rbind(c(0, 0, 0, 0, 0), c(1, 1, 1, 1, 1)))))
    
    data_temp <- MRData(covs = list(cum_mob = array(seq(0,1,length.out = 100))))
    
    model$attach_data(data_temp)
    
    
    model$fit_model()
    
    
    #Predict Spline
    
    global_pred <- expand.grid(location_id = unique(dt$location_id),
                               antigen = c("dtp", "mcv"),
                               cum_mob = c(unique(dt$cum_mob))) %>% as.data.table()
    
    dat_pred <- MRData()
    
    dat_pred$load_df(
      data = global_pred,
      col_covs=list("cum_mob")
    )
    
    global_pred$pred <- model$predict(dat_pred) %>% exp()
    
    global_pred <- merge(global_pred, hierarchy[,.(location_id, location_name)])
    global_pred <- merge(global_pred, 
                         dt[, .(location_id, cum_mob, obs = ratio, lower, upper, weight, antigen)],
                         by = c("cum_mob", "location_id"),
                         allow.cartesian = T, all.x = T)
    
    
    gg <- ggplot(global_pred, aes(x = cum_mob, y = pred))+
      geom_point(data = global_pred[!is.na(obs)], aes(y = obs, size = weight), alpha = 0.25, show.legend = F)+
      geom_line(color = "red")+
      labs(title = "Global splines: Both Antigens: disruption vs. mobility",
           x = "Cumulative Mobility",
           y = "Cumulative Vaccine Disruption Ratio")+
      scale_y_continuous(limits = c(0, 1.25), oob = squish)+
      theme_bw()
    
    print(gg)
    
    # cascading splines
    fit1 <- run_spline_cascade(
      stage1_model_object = model,
      df = dt,
      col_obs = "log_ratio",
      col_obs_se = "log_ratio_se",
      col_study_id = "loc_antigen",
      stage_id_vars = c("super_region_id","sr_antigen", "loc_antigen"),
      gaussian_prior = TRUE,
      inner_print_level = 0L,
      thetas = thetas_1,
      output_dir = work_dir,
      model_label = "stage1",
      overwrite_previous = TRUE
    )
    
    pred_dtp <- copy(mob[location_id %in% countries])
    pred_dtp[, antigen := "dtp"]
    
    pred_mcv <- copy(mob[location_id %in% countries])
    pred_mcv[, antigen := "mcv"]
    
    pred_dt <- rbind(pred_dtp, pred_mcv)
    
    # SR Preds
    sr_dt <- expand.grid(location_id = hierarchy[location_id == super_region_id, location_id],
                         cum_mob = seq(0,.7,0.01),
                         date = as.Date("2020-06-01"))
    
    sr_antigen_dt <- expand.grid(location_id = hierarchy[location_id == super_region_id, location_id],
                                 cum_mob = seq(0,.7,0.01),
                                 antigen = c("dtp", "mcv"),
                                 date = as.Date("2020-06-01"))
    
    ant_dt <- expand.grid(location_id = 1,
                          cum_mob = seq(0, .7, 0.01),
                          antigen = c("dtp", "mcv"),
                          date = as.Date("2020-06-01"))
    
    pred_dt <- rbindlist(list(pred_dt, sr_dt, sr_antigen_dt, ant_dt), fill = T)
    
    pred_dt <- merge(pred_dt, hierarchy[, .(location_id, super_region_id, super_region_name)], by = "location_id", all.x = T)
    
    pred_dt[, loc_antigen := paste0(location_id, "_", antigen)]
    pred_dt[, sr_antigen := paste0(super_region_id, "_", antigen)]
    
    preds1 <- predict_spline_cascade(
      fit = fit1,
      newdata = pred_dt
    ) %>% as.data.table()
    
    preds1[, pred:= exp(pred)]
    preds1[, location_id := as.integer(location_id)]
    
    gg <- ggplot(preds1[grep("r", cascade_prediction_id)], aes(x = cum_mob, y = pred))+
      geom_point(data = dt, aes(y = ratio, size = weight, color = antigen), alpha = 0.2, show.legend = F)+
      geom_line(aes(color = antigen, linetype = "Region"))+
      geom_line(data = global_pred, aes(linetype = "Global"))+
      facet_wrap(~super_region_name)+
      labs(title = "Regional Fits",
           x = "Cumulative Mobility Disruption",
           y = "Cumulative Vacccine Disruption",
           color = "")+
      scale_linetype_manual(values = c("Region" = "solid", "Global" = "dashed"))+
      scale_y_continuous(limits = c(0, 1.25), oob = squish)+
      theme_bw()
    
    if("antigen__dtp" %in% preds1$cascade_prediction_id){
      gg <- gg+geom_line(data = preds1[cascade_prediction_id %in% c("antigen__dtp", "antigen__mcv"), .(antigen, cum_mob, pred)],
                         aes(color = antigen, linetype = "Global"))
    }
    
    print(gg)
    
    preds1[, location_id := as.integer(location_id)]
    
    preds1 <- merge(preds1, hierarchy[, .(location_id, location_name)], by = "location_id")
    
    dt[,location_name := NULL]
    dt <- merge(dt, hierarchy[, .(location_id, location_name)], by = "location_id")
    
    for(r_name in unique(dt$super_region_name)){
      
      gg <- ggplot(preds1[grepl("loc", cascade_prediction_id) & super_region_name == r_name], 
                   aes(x = cum_mob, y = pred))+
        geom_point(data = dt[super_region_name == r_name], aes(y = ratio, size = weight, color = antigen), 
                   alpha = 0.4, show.legend = F)+
        geom_line(data = preds1[location_name == r_name, .(cum_mob, pred, antigen)],
                  aes(color = antigen, linetype = "Region"))+
        geom_line(aes(color = antigen, linetype = "Country"))+
        geom_line(data = unique(global_pred[, .(cum_mob, pred)]), aes(linetype = "Global"))+
        facet_wrap(~location_name, scales = "free")+
        labs(title = paste0("Country Fits: ", r_name),
             x = "Cumulative Mobility Disruption",
             y = "Cumulative Vacccine Disruption",
             color = "",
             linetype = "")+
        scale_linetype_manual(values = c("Country" = "solid", "Region" = "dashed", "Global" = "dotted"))+
        theme_bw()
      
      print(gg)
      
      
    }
    
    
    
    # Stage 2: Residuals vs. Time ---------------------------------------------
    
    # logit transform data
    
    resid_dt <- copy(dt)
    
    # For sources reporting weekly, only take 1 cumulative ratio per month, the date closest to the end of the month
    resid_dt[, nearest_month_end := last_days[which.min(abs(last_days-date))], by = 1:nrow(resid_dt)]
    
    resid_dt <- resid_dt[date %in% c(last_days, end_week_last_days)]
    
    resid_dt[, month := month(nearest_month_end)]
    
    # aggregate subnats into 1 national with population weights
    tmp_resid_dt <- resid_dt[grepl("_", ihme_loc_id)]
    tmp_resid_dt[, child_id := tstrsplit(ihme_loc_id, "_", keep = 2)]
    tmp_resid_dt[, child_id := as.integer(child_id)]
    
    tmp_resid_dt <- merge(tmp_resid_dt, pop[, .(location_id, parent_pop = population)], by = "location_id")
    tmp_resid_dt <- merge(tmp_resid_dt, pop[, .(location_id, child_pop = population)], by.x = "child_id", by.y = "location_id")
    
    tmp_resid_dt <- tmp_resid_dt[, .(ratio = weighted.mean(ratio, weights = child_pop), 
                                     weight = sum(weight),
                                     child_pop = sum(child_pop), 
                                     source = "Aggregated Sources"),
                                 by = c("location_id", "month", "super_region_id", "super_region_name", "location_name", 
                                        "parent_pop", "nearest_month_end", "antigen")]
    
    # only keep locations where we have subnats to account for at least half of the population (US and UK, but not Spain)
    tmp_resid_dt <- tmp_resid_dt[child_pop/parent_pop > 0.5]
    
    # recalculate log ratio and standard error from weighted mean of ratio and sum of inverse variance
    tmp_resid_dt[, log_ratio := log(ratio)]
    tmp_resid_dt[, log_ratio_se := 1/sqrt(weight)]
    
    tmp_resid_dt[, lower := exp(log_ratio - 1.96*log_ratio_se)]
    tmp_resid_dt[, upper := exp(log_ratio + 1.96*log_ratio_se)]
    
    resid_dt <- rbind(resid_dt[!grepl("_", ihme_loc_id)],
                      tmp_resid_dt, fill = T)
    
    
    # first convert to linear using MRBRT function
    resid_dt[, c("linear_ratio", "linear_ratio_se") := log_to_linear(array(log_ratio), array(log_ratio_se))]
    
    resid_dt[, date:= nearest_month_end]
    
    # merge on phase 1 model estimates
    resid_dt <- merge(resid_dt, preds1[, .(location_id, date, pred, antigen)], by = c("location_id", "date", "antigen"), all.x = T)
    
    # scale any values close to or greater than 1
    # if the ratio or the prediction is greater than 1-offset, choose a scalar to add to each value such that 
    #   neigher is greater than 1-offset. 
    
    resid_dt[, scalar:= min(c(0, (1-offset)-linear_ratio, (1-offset)-pred)), by = 1:nrow(resid_dt)]
    resid_dt[, linear_ratio := linear_ratio + scalar]
    resid_dt[, pred := pred + scalar]
    
    resid_dt[, logit_pred := logit(pred)]
    
    # convert to logit using MRBRT function
    resid_dt[, c("logit_ratio", "logit_ratio_se") := linear_to_logit(array(linear_ratio), array(linear_ratio_se))]
    
    #the SE values are extreme, so I set caps
    lower_bound <- quantile(resid_dt$logit_ratio_se, 0.25)
    upper_bound <- quantile(resid_dt$logit_ratio_se, 0.75)
    
    resid_dt[logit_ratio_se < lower_bound, logit_ratio_se := lower_bound]
    resid_dt[logit_ratio_se > upper_bound, logit_ratio_se := upper_bound]
    
    
    # Calculate residual
    resid_dt[, logit_resid := logit_ratio - logit_pred]
    
    # cap on residual to prevent extremes at 90th percentile
    resid_cap <- quantile(resid_dt$logit_resid, 0.90)
    resid_dt[logit_resid>resid_cap, logit_resid := resid_cap]
    resid_dt[logit_resid<(-1 * resid_cap), logit_resid := -1 * resid_cap]
    
    resid_dt[, loc_antigen := paste0(location_id, "_", antigen)]
    resid_dt[, sr_antigen := paste0(super_region_id, "_", antigen)]
    
    
    dat_2 <- MRData()
    dat_2$load_df(
      data = resid_dt,  col_obs = "logit_resid", col_obs_se = "logit_ratio_se",
      col_covs = list("month"), 
      col_study_id = "loc_antigen"
    )
    
    best_knots2 <- (quantile(resid_dt$month, c(0, 1/4, 2/4, 3/4, 1))-3)/9
    best_knots2[5] <- 1
    
    
    model2 <- MRBRT(data = dat_2,
                    cov_models = list(LinearCovModel("intercept", use_re = T),
                                      LinearCovModel("month",
                                                     use_spline = TRUE,
                                                     use_re = FALSE,
                                                     spline_degree = 2L,
                                                     spline_knots_type = "domain",
                                                     spline_knots = best_knots2,
                                                     spline_r_linear = TRUE,
                                                     spline_l_linear = TRUE)),
                    inlier_pct = 0.9)
    
    data_temp <- MRData(covs = list(month = array(seq(3,12, by = 0.25))))
    model2$attach_data(data_temp)
    
    
    model2$fit_model()
    
    
    #Predict Spline
    
    global_pred2 <- expand.grid(location_id = unique(resid_dt$location_id),
                                month = seq(0,12,by = 0.25),
                                antigen = c("dtp", "mcv")) %>% as.data.table()
    
    dat_pred <- MRData()
    
    dat_pred$load_df(
      data = global_pred2,
      col_covs=list("month")
    )
    
    global_pred2$pred <- model2$predict(dat_pred, predict_for_study = F)
    
    global_pred2 <- merge(global_pred2, hierarchy[,.(location_id, location_name)])
    global_pred2 <- merge(global_pred2, 
                          resid_dt[, .(location_id, logit_resid, logit_ratio_se, month)],
                          by = c("month", "location_id"),
                          allow.cartesian = T, all.x = T)
    
    
    gg <- ggplot(global_pred2, aes(x = month, y = pred))+
      geom_jitter(data = global_pred2[!is.na(logit_resid)], aes(y = logit_resid, size = logit_ratio_se^(-2)), 
                  width = 0.2, height = 0.05 ,alpha = 0.4, show.legend = F)+
      geom_line(show.legend = F, color = "red")+
      labs(title = "Global splines: residuals vs. time",
           x = "Month",
           y = "Residual (logit space)")+
      geom_hline(aes(yintercept = 0), linetype = "dotted")+
      theme_bw()
    
    print(gg)
    
    final_slope <- model2$cov_models[[2]]$spline$design_dmat(12, order = 1) %*% model2$beta_soln[1:4]
    
    model2$cov_models[[2]]$prior_spline_maxder_gaussian[, 4] <- c(final_slope, 0.2)
    model2$cov_models[[2]]$prior_spline_maxder_gaussian[2, 1:3] <- c(0.2, 0.2, 0.2)
    
    model2$inlier_pct <- 1
    
    # cascading splines
    fit2 <- run_spline_cascade(
      stage1_model_object = model2,
      df = resid_dt,
      col_obs = "logit_resid",
      col_obs_se = "logit_ratio_se",
      col_study_id = "location_id",
      stage_id_vars = c("super_region_id","sr_antigen", "loc_antigen"),
      gaussian_prior = TRUE,
      inner_print_level = 0L,
      thetas = thetas_2,
      output_dir = work_dir,
      model_label = "resid",
      overwrite_previous = TRUE
    )
    
    pred_dt2 <- expand.grid(location_id = c(countries, unique(hierarchy$super_region_id)), 
                            month = seq(0,12,0.25)) %>% as.data.table()
    
    pred_dt2_ant <- expand.grid(location_id = c(countries, unique(hierarchy$super_region_id)), 
                                month = seq(0,12,0.25), 
                                antigen = c("dtp", "mcv")) %>% as.data.table()
    
    pred_dt2 <- rbind(pred_dt2_ant, pred_dt2, fill = T)
    
    pred_dt2 <- merge(pred_dt2, hierarchy[, .(location_id, super_region_id, super_region_name, location_name)], by = "location_id", all.x = T)
    
    pred_dt2[, sr_antigen := paste0(super_region_id, "_", antigen)]
    pred_dt2[, loc_antigen := paste0(location_id, "_", antigen)]
    
    preds2 <- predict_spline_cascade(
      fit = fit2,
      newdata = pred_dt2
    ) %>% as.data.table()
    
    gg <- ggplot(preds2[grep("r", cascade_prediction_id)], aes(x = month, y = pred))+
      geom_point(data = resid_dt, aes(y = logit_resid, size = logit_ratio_se^(-2), color = antigen), alpha = 0.2, show.legend = F)+
      geom_line(aes(color = antigen, linetype = "Region"))+
      geom_line(data = global_pred2, aes(linetype = "Global"))+
      facet_wrap(~super_region_name)+
      labs(title = "Regional Fits: Residual",
           x = "Month",
           y = "Residuals in Logit Space",
           color = "",
           linetype = "")+
      scale_linetype_manual(values = c("Region" = "solid", "Global" = "dashed"))+
      theme_bw()
    
    if("antigen__dtp" %in% preds2$cascade_prediction_id){
      gg <- gg+geom_line(data = preds2[cascade_prediction_id %in% c("antigen__dtp", "antigen__mcv"), .(antigen, month, pred)],
                         aes(color = antigen, linetype = "Global"))
    }
    
    print(gg)
    
    
    for(r_name in unique(dt[!is.na(super_region_name), super_region_name])){
      
      gg <- ggplot(preds2[grepl("loc", cascade_prediction_id) & super_region_name == r_name], 
                   aes(x = month, y = pred))+
        geom_point(data = resid_dt[super_region_name == r_name], 
                   aes(y = logit_resid, size = logit_ratio_se^(-2), color = antigen), 
                   alpha = 0.4, show.legend = F)+
        geom_line(data = preds2[location_name == r_name, .(month, pred, antigen)],
                  aes(color = antigen, linetype = "Region"))+
        geom_line(aes(color = antigen, linetype = "Country"))+
        geom_line(data = unique(global_pred2[, .(month, pred)]), aes(linetype = "Global"))+
        facet_wrap(~location_name, scales = "free")+
        labs(title = paste0("Country Fits Residual: ", r_name),
             x = "Month",
             y = "Residual (logit space)",
             color = "", linetype = "")+
        scale_linetype_manual(values = c("Country" = "solid", "Region" = "dashed", "Global" = "dotted"))+
        theme_bw()
      
      print(gg)
      
      
    }
    
    dev.off()
    
    preds2[month %in% c(3:12), date := timeLastDayInMonth(paste0("2020-", month, "-01")) %>% as.Date()]
    preds2[,location_id := as.integer(location_id)]
    
    pred <- merge(preds1, preds2, by = c("date", "location_id", "super_region_id", "super_region_name", "location_name", "antigen"))
    
    pred[pred.x > (1-offset), pred.x := (1-offset)]
    pred[, logit_pred := logit(pred.x) + pred.y]
    pred[, ratio := inv_logit(logit_pred)]
    
    # convert to instantaneous
    pred <- pred[order(location_id, month, antigen)]
    pred[, month_ratio := ((month-2)*ratio)-((month-3)*data.table::shift(ratio)), by = c("location_id", "location_name", "antigen")]
    pred[month == 3, month_ratio := ratio]
    
    pred[, month_ratio.x := ((month - 2)*pred.x)-((month -3)*data.table::shift(pred.x)), by = c("location_id", "location_name", "antigen")]
    pred[month == 3, month_ratio.x := pred.x]
    
    write.csv(pred, file.path(work_dir, "predictions.csv"), row.names = F)
    
    # Convert back to daily space
    
    dt <- dt[order(source, location_id, end_date, antigen)]
    dt[, source_count := 1:.N, by = c("source", "ihme_loc_id", "antigen")]
    
    dt[, inst_ratio := (source_count*ratio)-((source_count - 1)*data.table::shift(ratio)), by = c("source", "ihme_loc_id", "antigen")]
    
    dt[source_count == 1, inst_ratio := ratio]
    
    # apply max standard error by source
    dt[, log_se_inst := max(log_ratio_se), by = c("ihme_loc_id", "source_id", "antigen")]
    dt[, log_ratio_inst := log(inst_ratio)]
    
    dt[, lower_inst := exp(log_ratio_inst - 1.96*log_se_inst)]
    dt[, upper_inst := exp(log_ratio_inst + 1.96*log_se_inst)]
    
    resid_dt <- resid_dt[order(source, location_id, end_date, antigen)]
    resid_dt[, source_count := 1:.N, by = c("source", "location_id", "antigen")]
    
    resid_dt[, inst_ratio := (source_count*ratio)-((source_count - 1)*data.table::shift(ratio)), by = c("source", "location_id", "antigen")]
    
    resid_dt[source_count == 1, inst_ratio := ratio]
    
    # apply max standard error by source
    resid_dt[, log_se_inst := max(log_ratio_se), by = c("location_id", "source_id", "antigen")]
    resid_dt[, log_ratio_inst := log(inst_ratio)]
    
    resid_dt[, lower_inst := exp(log_ratio_inst - 1.96*log_se_inst)]
    resid_dt[, upper_inst := exp(log_ratio_inst + 1.96*log_se_inst)]
    
    if(uncertainty){
      
      for(ant in c("dtp", "mcv")){
        
        ratio <- paste0("vaccine_", ant)
        
        draw_dir <- file.path(work_dir, ratio, "draws")
        
        if(file.exists(draw_dir)){
          unlink(draw_dir, recursive = TRUE)
        }
        
        dir.create(draw_dir, recursive = T)
        
        final_pred_dt <- mob_nat[date %in% last_days][order(location_id, date)]
        final_pred_dt <- merge(final_pred_dt, hierarchy[, .(location_id, location_name, super_region_id, super_region_name)])
        
        write.csv(final_pred_dt, file.path(work_dir, "pred_frame.csv"), row.names = F)
        
        write.csv(dt, file.path(work_dir, "input_data.csv"), row.names = F)
        write.csv(resid_dt, file.path(work_dir, "input_data_phase_2.csv"), row.names = F)
        
        for(loc_id in unique(dt$location_id)){
          
          script <- file.path("FILEPATH/02b_spline_cascade_draws.R")
          
          args <- paste(n_draws, loc_id, FALSE, work_dir, ant, offset, region_weight)
          mem <- "-l m_mem_free=5G"
          fthread <- "-l fthread=1"
          runtime <- "-l h_rt=00:05:00"
          archive <- "-l archive=TRUE"
          jname <- paste0("-N ", ant, "_loc_", loc_id)
          
          system(paste("qsub",jname,mem,fthread,runtime,archive,project,"-q all.q",sge.output.dir,rshell,
                       "-i FILEPATH/ihme_rstudio_3631.img","-s",script,args))
          
        }
        
        for(sr_id in hierarchy[!(location_id %in% dt$location_id) & level == 3, unique(super_region_id)]){
          
          script <- file.path("FILEPATH/02b_spline_cascade_draws.R")
          
          args <- paste(n_draws, sr_id, TRUE, work_dir, ant, offset, region_weight)
          mem <- "-l m_mem_free=5G"
          fthread <- "-l fthread=1"
          runtime <- "-l h_rt=00:05:00"
          archive <- "-l archive=TRUE"
          jname <- paste0("-N ", ant, "_sr_", sr_id)
          
          system(paste("qsub",jname,mem,fthread,runtime,archive,project,"-q all.q",sge.output.dir,rshell,
                       "-i FILEPATH/ihme_rstudio_3631.img","-s",script,args))
          
        }
        
        job_hold(ant)
        
        # Retry any that failed:
        
        for(loc_id in unique(dt$location_id)){
          
          if(paste0(loc_id, "_draws.csv") %in% list.files(draw_dir)){}else{
            script <- file.path("FILEPATH/02b_spline_cascade_draws.R")
            
            args <- paste(n_draws, loc_id, FALSE, work_dir, ant, offset, region_weight)
            mem <- "-l m_mem_free=100G"
            fthread <- "-l fthread=1"
            runtime <- "-l h_rt=00:05:00"
            archive <- "-l archive=TRUE"
            jname <- paste0("-N ", ant, "_loc_", loc_id)
            
            system(paste("qsub",jname,mem,fthread,runtime,archive,project,"-q all.q",sge.output.dir,rshell,
                         "-i FILEPATH/ihme_rstudio_3631.img","-s",script,args))
          }
        }
        
        for(sr_id in hierarchy[!(location_id %in% dt$location_id) & level == 3, unique(super_region_id)]){
          
          if(paste0(sr_id, "_draws.csv") %in% list.files(draw_dir)){}else{
            script <- file.path("FILEPATH/02b_spline_cascade_draws.R")
            
            args <- paste(n_draws, sr_id, TRUE, work_dir, ant, offset, region_weight)
            mem <- "-l m_mem_free=100G"
            fthread <- "-l fthread=1"
            runtime <- "-l h_rt=00:05:00"
            archive <- "-l archive=TRUE"
            jname <- paste0("-N ", ant, "_sr_", sr_id)
            
            system(paste("qsub}",jname,mem,fthread,runtime,archive,project,"-q all.q",sge.output.dir,rshell,
                         "-i FILEPATH/ihme_rstudio_3631.img","-s",script,args))
          }
        }
        
        job_hold(ant)
        
        final_pred_draws <- lapply(list.files(draw_dir, full.names = T), fread) %>% rbindlist()
        
        write.csv(final_pred_draws, file.path(file.path(work_dir, ratio, "all_draws.csv")), row.names = F)
        
        final_pred_draws_long <- melt.data.table(final_pred_draws, 
                                                 id.vars = c("location_id", "month", "location_name", "super_region_name", "super_region_id"),
                                                 measure.vars = patterns("cum_draw_", "resid_draw_"),
                                                 variable.name = "draw",
                                                 value.name = c("cum_ratio_1", "resid"))
        
        final_pred_draws_long[, cum_ratio_1:= as.numeric(cum_ratio_1)]
        final_pred_draws_long[, resid:= as.numeric(resid)]
        
        final_pred_draws_long[, draw:= as.integer(as.factor(draw))]
        
        final_pred_draws_long[cum_ratio_1 > (1-offset), cum_ratio_1 := (1-offset)]
        final_pred_draws_long[, cum_ratio_2 := inv_logit(logit(cum_ratio_1) + resid)]
        
        final_pred_draws_long <- final_pred_draws_long[order(draw, location_id, month)]
        
        final_pred_draws_long[, month_ratio_1 := ((month-2)*cum_ratio_1)-((month-3)*data.table::shift(cum_ratio_1)), by = c("location_id", "draw")]
        final_pred_draws_long[, month_ratio_2 := ((month-2)*cum_ratio_2)-((month-3)*data.table::shift(cum_ratio_2)), by = c("location_id", "draw")]
        
        final_pred_draws_long[month == 3, month_ratio_1 := cum_ratio_1]
        final_pred_draws_long[month == 3, month_ratio_2 := cum_ratio_2]
        
        # set any monthly ratios below zero to be zero
        final_pred_draws_long[month_ratio_1<offset, month_ratio_1 := offset]
        final_pred_draws_long[month_ratio_2<offset, month_ratio_2 := offset]
        
        final_pred_sum <- final_pred_draws_long[, .(cum_ratio_1_mean = mean(cum_ratio_1),
                                                    cum_ratio_1_lower = quantile(cum_ratio_1, p = 0.025),
                                                    cum_ratio_1_upper = quantile(cum_ratio_1, p = 0.975),
                                                    resid_mean = mean(resid),
                                                    resid_lower = quantile(resid, p = 0.025),
                                                    resid_upper = quantile(resid, p = 0.975),
                                                    cum_ratio_2_mean = mean(cum_ratio_2),
                                                    cum_ratio_2_lower = quantile(cum_ratio_2, p = 0.025),
                                                    cum_ratio_2_upper = quantile(cum_ratio_2, p = 0.975),
                                                    month_ratio_1_mean = mean(month_ratio_1),
                                                    month_ratio_1_lower = quantile(month_ratio_1, p = 0.025),
                                                    month_ratio_1_upper = quantile(month_ratio_1, p = 0.975),
                                                    month_ratio_2_mean = mean(month_ratio_2),
                                                    month_ratio_2_lower = quantile(month_ratio_2, p = 0.025),
                                                    month_ratio_2_upper = quantile(month_ratio_2, p = 0.975)),
                                                by = c("location_id", "month", "location_name", "super_region_name", "super_region_id")] %>% unique()
        
        write.csv(final_pred_sum, file.path(file.path(work_dir, ratio, "summary_results.csv")), row.names = F)
      }
    }
  
    print(paste("Finished with iteration", i, "out of 5, month ", n, " out of 10:", Sys.time() - start_time))
  }
}




print(Sys.time() - start_time)





