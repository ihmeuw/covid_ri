#-------------------Header------------------------------------------------
# Author: NAME
# Date: 2/8/21
# Purpose: Read in coverage data and calculate coverage/disruptions from draws
#          Create Global and Super-Region aggregates
#          
# source("FILEPATH/summarize_results.R", echo=T)
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

packages <- c("data.table","magrittr","ggplot2", "RNetCDF")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Directories -------------------------------------------------------------

dir <- "VERSION_NAME"


input_dir <- paste0(j_root, "/Project/goalkeepers_2020/covid_interruption_covariates/cumulative_", dir)
output_dir <- paste0(j_root, "/temp/kcausey/vaccine_cumulative/", dir)

dir.create(output_dir, recursive = T)

source(file.path(central_lib, "FILEPATH/get_population.R"))
source(file.path(central_lib, "FILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)

out_data_path <- file.path(output_dir, "summary_results.csv")

offset <- 1e-3

logit <- function(x){log(x/(1-x))}
inv_logit <- function(x){exp(x)/(1+exp(x))}

# Pull disruption ---------------------------------------------------------

model_names <- c("vaccine_dtp", "vaccine_mcv")

model_results <- data.table()

for(model_name in model_names){
  tmp_results <- fread(paste0(input_dir, "/",model_name, "/all_draws.csv"))
  tmp_results[, model := model_name]
  model_results <- rbind(model_results, tmp_results)
}

# make long by draw
model_results_long <- melt.data.table(model_results, 
                                         id.vars = c("location_id", "month", "location_name", "super_region_name", "super_region_id","model"),
                                         measure.vars = patterns("cum_draw_", "resid_draw_"),
                                         variable.name = "draw",
                                         value.name = c("cum_ratio_1", "resid"))

model_results_long[, draw := as.integer(as.factor(draw))]
model_results_long[, cum_ratio_1 := as.numeric(cum_ratio_1)]
model_results_long[, resid := as.numeric(resid)]

model_results_long[cum_ratio_1 > (1-offset), cum_ratio_1 := (1-offset)]
model_results_long[, cum_ratio_2 := inv_logit(logit(cum_ratio_1) + resid)]

model_results_long <- model_results_long[order(model, draw, location_id, month)]

model_results_long[, month_ratio_1 := ((month-2)*cum_ratio_1)-((month-3)*data.table::shift(cum_ratio_1)), by = c("location_id", "draw", "model")]
model_results_long[, month_ratio_2 := ((month-2)*cum_ratio_2)-((month-3)*data.table::shift(cum_ratio_2)), by = c("location_id", "draw", "model")]

model_results_long[month == 3, month_ratio_1 := cum_ratio_1]
model_results_long[month == 3, month_ratio_2 := cum_ratio_2]

# add in values for Jan and Feb
tmp_dt <- copy(model_results_long[, .(location_id, location_name, super_region_name, super_region_id, model, draw, month_ratio_2 = 1)] %>% unique)

model_results_long <- rbind(model_results_long, tmp_dt, fill = T)
model_results_long[is.na(month), month := 1]
model_results_long <- rbind(model_results_long, tmp_dt, fill = T)
model_results_long[is.na(month), month := 2]

# set any monthly ratios below zero to be zero
model_results_long[month_ratio_1<offset, month_ratio_1 := offset]
model_results_long[month_ratio_2<offset, month_ratio_2 := offset]


# Pull Reference Coverage -----------------------------------------------------------
files <- list.files("FILEPATH", full.names = T)

dtp3 <- lapply(files, fread) %>% rbindlist()
dtp3 <- dtp3[year_id == 2020 & location_id %in% model_results$location_id & year_id == 2020]
dtp3 <- melt.data.table(dtp3, id.vars = c("location_id"), measure.vars = patterns("draw_"), value.name = "reference_coverage")
dtp3[, draw := as.integer(as.factor(variable))]
dtp3[, variable := NULL]
dtp3[, model := "vaccine_dtp"]



files <- list.files("FILEPATH", full.names = T)

mcv1 <- lapply(files, fread) %>% rbindlist()
mcv1 <- mcv1[year_id == 2020 & location_id %in% model_results$location_id & year_id == 2020]
mcv1 <- melt.data.table(mcv1, id.vars = c("location_id"), measure.vars = patterns("draw_"), value.name = "reference_coverage")
mcv1[, draw := as.integer(as.factor(variable))]
mcv1[, variable := NULL]
mcv1[, model := "vaccine_mcv"]

coverage <- rbind(mcv1, dtp3)


# Calculate disruption coverage -------------------------------------------

setnames(model_results_long, "month_ratio_2", "value")
model_results_long[, draw:= as.integer(draw)]

out_dt <- merge(model_results_long, coverage, by = c("location_id", "model", "draw"), all.x = T)

out_dt[, coverage := reference_coverage * value]

# Aggregate ---------------------------------------------------------------

# pull 2020 populations for 0-11 and 12-23 months
pop_dt <- get_population(age_group_id = c(28, 238), year_id = 2020, location_id = unique(out_dt$location_id), gbd_round_id = 7, decomp_step = "iterative", single_year_age = T)

pop_dt <- dcast.data.table(pop_dt, location_id ~ age_group_id, value.var = "population")

setnames(pop_dt, c("28", "238"), c("pop_0", "pop_1"))

# for MCV1 some countries vaccinate after 1st bday and others before. Decide which populations to use
mcv1_schedule <- readRDS("FILEPATH/vaccine_target.rds")[me_name == "vacc_mcv1" & age_cohort > 1]
mcv1_schedule <- merge(mcv1_schedule, hierarchy[,.(location_id, ihme_loc_id)], all.x = T)

# merge on populations 

agg_dt <- merge(out_dt, pop_dt, by = c("location_id"), all.x = T)

agg_dt[, population := pop_0]
agg_dt[model == "vaccine_mcv" & location_id %in% mcv1_schedule$location_id, population := pop_1]
agg_dt[, c("pop_0", "pop_1") := NULL]
agg_dt[, level := 2]

regions <- agg_dt[, .(value = weighted.mean(value, w = population),
                      coverage = weighted.mean(coverage, w = population),
                      reference_coverage = weighted.mean(reference_coverage, w = population),
                      population = sum(population), level = 1), by = c("super_region_id", "super_region_name", "model", "draw", "month")]

setnames(regions, c("super_region_id", "super_region_name"), c("location_id", "location_name"))
regions <- regions[!is.na(location_id)]


global <- agg_dt[, .(value = weighted.mean(value, w = population),
                     coverage = weighted.mean(coverage, w = population),
                     reference_coverage = weighted.mean(reference_coverage, w = population),
                     population = sum(population), level = 0), by = c("model", "draw", "month")]

global[, location_id := 1]
global[, location_name := "Global"]

out_dt <- rbindlist(list(agg_dt[, names(regions), with = F], regions, global), use.names = TRUE)


# Calculate missed doses
out_dt[, n_missed := (population)/12 * (1-coverage)]
out_dt[, n_missed_ref := (population)/12 * (1-reference_coverage)]

out_dt[, period := month.abb[month]]

out_dt_early<- out_dt[month %in% 3:6, 
                                .(value = mean(value), coverage = mean(coverage), reference_coverage = mean(reference_coverage), 
                                  period = "Mar - Jun", n_missed = sum(n_missed), n_missed_ref = sum(n_missed_ref)),
                                by = c("location_id", "model", "draw", "location_name", "level", "population")]

out_dt_late <- out_dt[month %in% 7:12, 
                            .(value = mean(value), coverage = mean(coverage), reference_coverage = mean(reference_coverage), 
                              period = "Jul - Dec", n_missed = sum(n_missed), n_missed_ref = sum(n_missed_ref)),
                            by = c("location_id", "model", "draw", "location_name", "level", "population")]

out_dt_full <- out_dt[month %in% 1:12, 
                      .(value = mean(value), coverage = mean(coverage), reference_coverage = mean(reference_coverage), 
                        period = "Jan - Dec", n_missed = sum(n_missed), n_missed_ref = sum(n_missed_ref)),
                      by = c("location_id", "model", "draw", "location_name",  "level", "population")]

out_dt[, month := NULL]

out_dt <- rbindlist(list(out_dt, out_dt_early, out_dt_late, out_dt_full), use.names = T)

out_dt[, n_missed_covid := n_missed - n_missed_ref]

write.csv(out_dt, file.path(output_dir, "draws_results_agg.csv"), row.names = F)


# Summaries from draws ----------------------------------------------------

out_dt_sum <- out_dt[, .(mean = mean(value),
                         lower = quantile(value, 0.025),
                         upper = quantile(value, 0.975),
                         mean_coverage = mean(coverage),
                         lower_coverage = quantile(coverage, 0.025),
                         upper_coverage = quantile(coverage, 0.975),
                         mean_coverage_ref = mean(reference_coverage),
                         lower_coverage_ref = quantile(reference_coverage, 0.025),
                         upper_coverage_ref = quantile(reference_coverage, 0.975),
                         mean_n_missed = mean(n_missed),
                         lower_n_missed = quantile(n_missed, 0.025),
                         upper_n_missed = quantile(n_missed, 0.975),
                         mean_n_missed_ref = mean(n_missed_ref),
                         lower_n_missed_ref = quantile(n_missed_ref, 0.025),
                         upper_n_missed_ref = quantile(n_missed_ref, 0.975),
                         mean_n_missed_covid = mean(n_missed_covid),
                         lower_n_missed_covid = quantile(n_missed_covid, 0.025),
                         upper_n_missed_covid = quantile(n_missed_covid, 0.975)),
                     by = c("location_id", "model", "location_name", "level", "period")]


write.csv(out_dt_sum, file.path(output_dir, "summary_results_agg.csv"), row.names = F)


# Reshape and output results ----------------------------------------------

out_dt <- fread(file.path(output_dir, "summary_results_agg.csv"))

out_dt[, model := gsub("vaccine_", "", model)]

format_dt <- dcast.data.table(out_dt[period == "Jan - Dec"], location_id + location_name + level ~ model,
                              value.var = c("mean", "lower", "upper", 
                                            "mean_coverage", "lower_coverage", "upper_coverage", 
                                            "mean_coverage_ref", "lower_coverage_ref", "upper_coverage_ref",
                                            "mean_n_missed", "lower_n_missed", "upper_n_missed",
                                            "mean_n_missed_ref", "lower_n_missed_ref", "upper_n_missed_ref",
                                            "mean_n_missed_covid", "lower_n_missed_covid", "upper_n_missed_covid"))

n <- 1

format_dt <- format_dt[, .(location_name, location_id, level,
                           dtp_ref_cov = 
                             
                             paste0((mean_coverage_ref_dtp * 100) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_coverage_ref_dtp * 100) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_coverage_ref_dtp * 100) %>% round(digits = n) %>% format(nsmall = n), ")"),
                          
                           dtp_disrupt =

                             paste0(((1-mean_dtp) * 100) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    ((1-upper_dtp) * 100) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    ((1-lower_dtp) * 100) %>% round(digits = n) %>% format(nsmall = n), ")"),

                           dtp_pand_cov =

                             paste0((mean_coverage_dtp * 100) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_coverage_dtp * 100) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_coverage_dtp * 100) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           dtp_n_missed_ref = 
                             
                             paste0((mean_n_missed_ref_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_n_missed_ref_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_n_missed_ref_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           dtp_n_missed_covid = 
                             
                             paste0((mean_n_missed_covid_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_n_missed_covid_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_n_missed_covid_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           dtp_n_missed = 
                             
                             paste0((mean_n_missed_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_n_missed_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_n_missed_dtp/1e6) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           mcv_ref_cov = 
                             
                             paste0((mean_coverage_ref_mcv * 100) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_coverage_ref_mcv * 100) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_coverage_ref_mcv * 100) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           mcv_disrupt =
                             
                             paste0(((1-mean_mcv) * 100) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    ((1-upper_mcv) * 100) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    ((1-lower_mcv) * 100) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           mcv_pand_cov =
                             
                             paste0((mean_coverage_mcv * 100) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_coverage_mcv * 100) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_coverage_mcv * 100) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           mcv_n_missed_ref = 
                             
                             paste0((mean_n_missed_ref_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_n_missed_ref_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_n_missed_ref_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           mcv_n_missed_covid = 
                             
                             paste0((mean_n_missed_covid_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_n_missed_covid_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_n_missed_covid_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), ")"),
                           
                           mcv_n_missed = 
                             
                             paste0((mean_n_missed_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), "\n(",
                                    (lower_n_missed_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), ", ",
                                    (upper_n_missed_mcv/1e6) %>% round(digits = n) %>% format(nsmall = n), ")"))]

format_dt <- format_dt[order(level, location_name)]

write.csv(format_dt, out_data_path, row.names = F)



