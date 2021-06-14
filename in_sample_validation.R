
#-------------------Header------------------------------------------------
# Author: NAME
# Date: 3/8/2021
# Purpose: Create in-sample validation analysis and plot
#          
# source("FILEPATH/in_sample_validation.R", echo=T)
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

color_scheme <- c("DTP3" = "#a8216b",
                  "MCV1" = "#57c1e4")

work_dir <- file.path(j_root, "FILEPATH", paste0("cumulative_", model_name))
output_dir <- file.path(j_root, "FILEPATH", model_name)
dir.create(output_dir)

source(file.path(central_lib, "FILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)

antigens <- c("dtp", "mcv")

pred_dt1 <- data.table()
pred_dt2 <- data.table()



for(ant in antigens){
  
  dt1 <- fread(file.path(work_dir, "stage1", "input_data.csv"))[antigen == ant]
  dt2 <- fread(file.path(work_dir, "resid", "input_data.csv"))[antigen == ant]
  
  dt1[, model := ant]
  dt2[, model := ant]
  
  pred <- fread(paste0(work_dir, "/vaccine_", ant, "/summary_results.csv"))
  
  dt1 <- merge(dt1, pred[, .(location_id, month, pred_1 = cum_ratio_1_mean)], by = c("location_id", "month"))
  dt2 <- merge(dt2, pred[, .(location_id, month, pred_2 = cum_ratio_2_mean)], by = c("location_id", "month"))
  
  pred_dt1 <- rbind(pred_dt1, dt1)
  pred_dt2 <- rbind(pred_dt2, dt2)
}

pred_dt1[, date := as.Date(date)]

pred_dt1[, max_date := max(date), by = c("model", "ihme_loc_id", "source")]

pred_dt1 <- pred_dt1[, .(location_id, date, source_id, location_name, ihme_loc_id, source,
                         ratio, log_ratio_se, weight, lower, upper, cum_mob, month, super_region_name,
                         pred = pred_1, model, phase = "Step 1", max_date)]

pred_dt2[, date := as.Date(date)]

pred_dt2[, max_date := max(date), by = c("model", "location_id", "source")]

pred_dt2 <- pred_dt2[, .(location_id, date, source_id, location_name, ihme_loc_id, source,
                         ratio, log_ratio_se, weight, lower, upper, cum_mob, month, super_region_name,
                         pred = pred_2, model, phase = "Step 2", max_date)]

dt <- rbind(pred_dt1, pred_dt2)

dt <- dt[date == max_date]
dt[, max_date := NULL]

dt[model == "dtp", antigen := "DTP3"]
dt[model == "mcv", antigen := "MCV1"]


dt[, resid := pred - ratio]
dt[, abs_resid := abs(pred - ratio)]

dt <- dt[!grepl("_", ihme_loc_id)]


model_fit <- dt[, .(mean_error = mean(resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    weighted_mean_error = weighted.mean(resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    mean_absolute_error = mean(abs_resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    weighted_mean_absolute_error = weighted.mean(abs_resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                    RMSE = RMSE(pred, ratio) %>% signif(3) %>% sprintf(fmt = "%.4f")), by = c("antigen", "phase")]

model_fit

model_fit_pos <- dt[ratio<1, .(mean_error = mean(resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               weighted_mean_error = weighted.mean(resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               mean_absolute_error = mean(abs_resid) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               weighted_mean_absolute_error = weighted.mean(abs_resid, weight) %>% signif(3) %>% sprintf(fmt = "%.4f"),
                               RMSE = RMSE(pred, ratio) %>% signif(3) %>% sprintf(fmt = "%.4f")), by = c("antigen", "phase")]


gg <- ggplot(dt, aes(x = ratio, y = pred))+
  geom_abline(aes(slope = 1, intercept = 0), color = "grey50")+
  geom_point(aes(size = weight, color = antigen),alpha = 0.2)+
  geom_text(data = model_fit, aes(x = 1, y = 0.4, 
                                  label = paste0(paste0("wME = ", weighted_mean_error, 
                                                        "\n","wMAE = ", weighted_mean_absolute_error, 
                                                        "\n","RMSE = ", RMSE))))+
  facet_grid(phase~antigen)+
  scale_size_continuous(breaks = c(10, 100, 1000, 10000, 100000))+
  labs(y = "Modeled Estimate",
       x = "Data",
       size = "Weight")+
  coord_fixed(xlim = c(0.30,1.35), ylim = c(0.30, 1.35), expand = 0)+
  scale_color_manual(values = color_scheme, guide = F)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))

pdf(file.path(output_dir, "in_sample_validation_plot.pdf"), width = 10, height = 6)

print(gg)

dev.off()





