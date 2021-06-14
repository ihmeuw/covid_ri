
#-------------------Header------------------------------------------------
# Author: NAME
# Date: 3/8/2021
# Purpose: Create out of sample validation analysis and plot for paper
#          
# source("FILEPATH/out_of_sample_validation.R", echo=T)
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

packages <- c("data.table","magrittr","DescTools", "ggplot2", "ggrepel", "spatstat")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Directories -------------------------------------------------------------

library(mrbrt001, lib.loc = "FILEPATH")

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

model_name <- "VERSION_NAME"

compare_version <- "VERSION_NAME"

color_scheme <- c("DTP3" = "#a8216b",
                  "MCV1" = "#57c1e4")

work_dir <- file.path(j_root, "FILEPATH", paste0("cumulative_", model_name))
output_dir <- work_dir

compare_dir <- file.path(j_root, "FILEPATH", paste0("cumulative_", compare_version))

source(file.path(central_lib, "FILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)

color_scheme_6 <- c("#a8216b",
                    "#57c1e4",
                    "#ec1b4b",
                    "#2e9598",
                    "#f7db69",
                    "#f26a44",
                    "#0A2F51")

antigens <- c("dtp", "mcv")

dt <- fread(file.path(work_dir, "input_data.csv"))

holdout_locs <- fread(file.path(work_dir, "holdout_locs.csv"))

dt <- merge(dt, holdout_locs[, .(location_id, holdout)], by = "location_id")

# create a file of predictions
pred_dt <- data.table()

for(i in 1:5){
  for(n in c(1:10)){
    for(ant in antigens){
      pred_tmp <- fread(paste0(work_dir, "/months_", n, "/holdout_", i, "/vaccine_", ant, "/summary_results.csv"))
      
      pred_tmp[, holdout := i]
      pred_tmp[, n_months := n]
      pred_tmp[, antigen := ant]
      
      pred_dt <- rbind(pred_dt, pred_tmp)
    }
  }
}

dtp <- dt[dtp == 1]
mcv <- dt[mcv == 1]

dtp[, antigen := "dtp"]
mcv[, antigen := "mcv"]

dt <- rbind(dtp, mcv)

dt <- merge(dt, pred_dt[, .(location_id, month, antigen, holdout, n_months, pred = cum_ratio_2_mean)], 
            by = c("location_id", "month", "antigen", "holdout"),
            allow.cartesian = T)

dt[, date := as.Date(date)]

dt[, max_date := max(date), by = c("antigen", "ihme_loc_id", "source")]

dt <- dt[date == max_date]
dt[, max_date := NULL]

# only keep locations with data remaining
dt[, max_month := max(month), by = c("antigen", "ihme_loc_id", "source")]
dt <- dt[n_months < max_month - 2]

dt[antigen == "dtp", vacc := "DTP3"]
dt[antigen == "mcv", vacc := "MCV1"]


dt[, resid := pred - ratio]
dt[, abs_resid := abs(pred - ratio)]

dt <- dt[!grepl("_", ihme_loc_id)]


model_fit <- dt[, .(mean_error_c = mean(resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    weighted_mean_error_c = weighted.mean(resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    mean_absolute_error_c = mean(abs_resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    weighted_mean_absolute_error_c = weighted.mean(abs_resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    RMSE_c = RMSE(pred, ratio) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    mean_error = mean(resid),
                    weighted_mean_error = weighted.mean(resid, weight),
                    mean_absolute_error = mean(abs_resid),
                    weighted_mean_absolute_error = weighted.mean(abs_resid, weight),
                    RMSE = RMSE(pred, ratio),
                    n = .N), 
                by = c("vacc", "n_months")]

model_fit

model_fit_pos <- dt[ratio<1, .(mean_error_c = mean(resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               weighted_mean_error_c = weighted.mean(resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               mean_absolute_error_c = mean(abs_resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               weighted_mean_absolute_error_c = weighted.mean(abs_resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               RMSE_c = RMSE(pred, ratio) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               mean_error = mean(resid),
                               weighted_mean_error = weighted.mean(resid, weight),
                               mean_absolute_error = mean(abs_resid),
                               weighted_mean_absolute_error = weighted.mean(abs_resid, weight),
                               RMSE = RMSE(pred, ratio),
                               n = .N), 
                    by = c("vacc", "n_months")]

model_fit_long <- melt(model_fit, id.vars = c("vacc", "n_months"),
                       measure.vars = c("n", "weighted_mean_error", "weighted_mean_absolute_error", "RMSE"))

model_fit_long[variable == "weighted_mean_error", label := "wME"]
model_fit_long[variable == "weighted_mean_absolute_error", label := "wMAE"]
model_fit_long[variable == "RMSE", label := "RMSE"]
model_fit_long[variable == "n", label := "N"]

model_fit_long[, label := factor(label, levels = c("N", "ME", "wME", "MAE", "wMAE", "RMSE"))]

model_fit_pos_long <- melt(model_fit_pos, id.vars = c("vacc", "n_months"),
                       measure.vars = c("n", "weighted_mean_error", "weighted_mean_absolute_error", "RMSE"))

model_fit_pos_long[variable == "weighted_mean_error", label := "wME"]
model_fit_pos_long[variable == "weighted_mean_absolute_error", label := "wMAE"]
model_fit_pos_long[variable == "RMSE", label := "RMSE"]
model_fit_pos_long[variable == "n", label := "N"]

model_fit_pos_long[, label := factor(label, levels = c("N", "ME", "wME", "MAE", "wMAE", "RMSE"))]

# Plot model fit stats

gg <- ggplot(model_fit_long, aes(x = n_months, y = value, color = variable))+
  geom_point(show.legend = F)+
  geom_point(data = model_fit_pos_long, alpha = 0, show.legend = F)+
  geom_line(show.legend = F)+
  geom_hline(aes(yintercept = 0), linetype = "dotted")+
  facet_grid(label~vacc, scales = "free")+
  scale_x_continuous(breaks = c(1:9))+
  scale_color_manual(values = color_scheme_6,
                     breaks = c("n", "mean_error", "weighted_mean_error", "mean_absolute_error", "weighted_mean_absolute_error", "RMSE"),
                     labels = c("N", "ME", "wME", "MAE", "wMAE", "RMSE"))+
  labs(title = "Out-of-sample validation: All CDRs",
       y = "",
       x = "Number of months witheld",
       color = "Metric")+
  theme_bw()

print(gg)


gg <- ggplot(model_fit_pos_long, aes(x = n_months, y = value, color = variable))+
  geom_point(show.legend = F)+
  geom_point(data = model_fit_long, alpha = 0, show.legend = F)+
  geom_line(show.legend = F)+
  geom_hline(aes(yintercept = 0), linetype = "dotted")+
  facet_grid(label~vacc, scales = "free")+
  scale_x_continuous(breaks = c(1:9))+
  scale_color_manual(values = color_scheme_6,
                     breaks = c("n", "weighted_mean_error", "weighted_mean_absolute_error", "RMSE"),
                     labels = c("N", "wME", "wMAE", "RMSE"))+
  labs(title = "Out-of-sample validation: CDRs < 1",
       y = "",
       x = "Number of months witheld",
       color = "Metric")+
  theme_bw()

print(gg)



# Plot data against predictions by months of holdout

gg <- ggplot(dt[vacc == "DTP3"], aes(x = ratio, y = pred))+
  geom_abline(aes(slope = 1, intercept = 0), color = "grey50")+
  geom_point(aes(size = weight, color = vacc),alpha = 0.2, show.legend = F)+
  geom_text(data = model_fit[vacc == "DTP3"], aes(x = 1, y = 0.5, 
                                  label = paste0("wME = ", weighted_mean_error_c, 
                                                 "\n","wMAE = ", weighted_mean_absolute_error_c, 
                                                 "\n","RMSE = ", RMSE_c,
                                                 "\n","n = ", n)),
            size = 3)+
  facet_wrap(~n_months)+
  labs(title = "DTP3",
       y = "Modeled Estimate",
       x = "Data")+
  coord_fixed(xlim = c(0.30,1.35), ylim = c(0.30, 1.35), expand = 0)+
  scale_color_manual(values = color_scheme, guide = F)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))

print(gg)

gg <- ggplot(dt[vacc == "MCV1"], aes(x = ratio, y = pred))+
  geom_abline(aes(slope = 1, intercept = 0), color = "grey50")+
  geom_point(aes(size = weight, color = vacc),alpha = 0.2, show.legend = F)+
  geom_text(data = model_fit[vacc == "MCV1"], aes(x = 1, y = 0.5, 
                                  label = paste0("wME = ", weighted_mean_error_c, 
                                                 "\n","wMAE = ", weighted_mean_absolute_error_c, 
                                                 "\n","RMSE = ", RMSE_c,
                                                 "\n","n = ", n)),
            size = 3)+
  facet_wrap(~n_months)+
  labs(title = "MCV1",
       y = "Modeled Estimate",
       x = "Data")+
  coord_fixed(xlim = c(0.30,1.35), ylim = c(0.30, 1.35), expand = 0)+
  scale_color_manual(values = color_scheme, guide = F)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))

print(gg)



# Analysis holding out all months of data-------------------------------------------

dt <- rbind(dtp, mcv)

dt <- merge(dt, pred_dt[, .(location_id, month, antigen, holdout, n_months, 
                            pred = cum_ratio_2_mean, pred_lower = cum_ratio_2_lower, pred_upper = cum_ratio_2_upper)], 
            by = c("location_id", "month", "antigen", "holdout"),
            allow.cartesian = T)

dt[, date := as.Date(date)]

dt[, max_date := max(date), by = c("antigen", "ihme_loc_id", "source")]

dt <- dt[date == max_date]
dt[, max_date := NULL]

dt[antigen == "dtp", vacc := "DTP3"]
dt[antigen == "mcv", vacc := "MCV1"]

dt <- dt[n_months == 10]

dt[, resid := pred - ratio]
dt[, abs_resid := abs(pred - ratio)]

dt <- dt[!grepl("_", ihme_loc_id)]


model_fit <- dt[, .(mean_error_c = mean(resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    weighted_mean_error_c = weighted.mean(resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    mean_absolute_error_c = mean(abs_resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    weighted_mean_absolute_error_c = weighted.mean(abs_resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    RMSE_c = RMSE(pred, ratio) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    mean_error = mean(resid),
                    weighted_mean_error = weighted.mean(resid, weight),
                    mean_absolute_error = mean(abs_resid),
                    weighted_mean_absolute_error = weighted.mean(abs_resid, weight),
                    RMSE = RMSE(pred, ratio),
                    n = .N), 
                by = c("vacc", "n_months")]

model_fit

model_fit_pos <- dt[ratio<1, .(mean_error = mean(resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               weighted_mean_error = weighted.mean(resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               mean_absolute_error = mean(abs_resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               weighted_mean_absolute_error = weighted.mean(abs_resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               RMSE = RMSE(pred, ratio) %>% signif(3) %>% sprintf(fmt = "%.4f")), by = c("vacc", "n_months")]

model_fit_long <- melt(model_fit, id.vars = c("vacc", "n_months", "n"),
                       measure.vars = c("mean_error", "weighted_mean_error", "mean_absolute_error", "weighted_mean_absolute_error", "RMSE"))


# Plot data against predictions dropping all months

gg <- ggplot(dt, aes(x = ratio, y = pred))+
  geom_abline(aes(slope = 1, intercept = 0), color = "grey50")+
  geom_point(aes(size = weight, color = vacc),alpha = 0.2, show.legend = F)+
  geom_text(data = model_fit, aes(x = 1, y = 0.5, 
                                                  label = paste0("wME = ", weighted_mean_error_c, 
                                                                 "\n","wMAE = ", weighted_mean_absolute_error_c, 
                                                                 "\n","RMSE = ", RMSE_c)),
            size = 3)+
  facet_wrap(~vacc)+
  labs(y = "Modeled estimate",
       x = "Data")+
  coord_fixed(xlim = c(0.30,1.35), ylim = c(0.30, 1.35), expand = 0)+
  scale_color_manual(values = color_scheme, guide = F)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))

print(gg)











