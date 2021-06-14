
#-------------------Header------------------------------------------------
# Author: NAME
# Date: 12/4/20
# Purpose: Read in vaccine extractions and reformat for MRBRT, cumulative space
#          
# source("FILEPATH", echo=T)
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

packages <- c("data.table","magrittr","ggplot2", "lubridate", "gridExtra", "grid", "stringr", 
              "msm", "wCorr", "timeDate", "dplyr")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}


# Directories -------------------------------------------------------------

data_dir <- paste0(Sys.Date(), "_drop_last_month") # Sys.Date() #

work_dir <- file.path(j_root, "FILEPATH")

admin_data_path <- file.path(work_dir, "COVID_immunization_data_available_by_country_20210404.csv") # update if admin data changes
who_data_path <- file.path(work_dir, "Copy_of_WHOAdminData_FullDataset_2020-11-13.csv")
paho_data_path <- file.path(work_dir, "paho_3_14.csv")
afro_data_path <- file.path(work_dir, "WHO_AFRO_3_24.csv")
emro_data_path <- file.path(work_dir, "emro_3_25.csv")
jrf_data_path <- file.path(work_dir, "jrf_2019_annual_doses.csv")

out_dir <- file.path(work_dir, "cumulative", data_dir)
dir.create(out_dir, recursive = T)

source(file.path(central_lib, "FILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)

last_days <- as.Date(timeLastDayInMonth(paste0("2020-", 1:12, "-01")))

end_week_last_days <-  as.Date(c("2020-02-01", "2020-02-29", "2020-03-28", "2020-05-02", "2020-05-30", "2020-06-27", 
                         "2020-08-01", "2020-08-29", "2020-10-03", "2020-10-31", "2020-11-28", "2021-01-02"))
# Read Data --------------------------------------------------------------

dt <- fread(admin_data_path)

dt <- merge(dt, hierarchy[, .(ihme_loc_id, location_id)], all.x=T)

who <- fread(who_data_path)
paho <- fread(paho_data_path)
afro <- fread(afro_data_path)
emro <- fread(emro_data_path)

jrf <- fread(jrf_data_path)

# drop problematic ghana points
afro <- afro[!(ISO3_Code == "GHA") ]

who <- rbind(paho, who[!(ISO3_Code %in% paho$ISO3_Code)], fill = T, use.names = T)

# pull afro from new when available and old data when not
who <- rbind(afro, who, use.names = T, fill = T)
who[, count:= 1:.N, by = c("ISO3_Code", "Month", "Year", "Antigen")]
who <- who[count == 1]

# pull emro from new when available and old data when not
who <- rbind(emro, who, use.names = T, fill = T)
who[, count:= 1:.N, by = c("ISO3_Code", "Month", "Year", "Antigen")]
who <- who[count == 1]

who <- merge(who, hierarchy[,.(ISO3_Code = ihme_loc_id, ihme_loc_id, location_id, location_name)], all.x=T)

# WHO ---------------------------------------------------------------

who <- who[!is.na(location_id) & Antigen %in% c("DTP3", "MCV1")]

who[, n_doses := gsub(" ", "", Number_of_doses_administered)]
who[, n_doses := gsub(",", "", n_doses)]
who[, n_doses := as.integer(n_doses)]

who <- who[, .(location_id, location_name, ihme_loc_id, year_id = Year, month_start = tolower(substring(Month, 1, 3)), antigen = tolower(Antigen), n_doses)] %>% unique()

who_dt <- dcast.data.table(who, ... ~ year_id, value.var = "n_doses")
setnames(who_dt, c("2019", "2020"), c("reference", "pandemic"))

who_dt[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(who_dt)]


# Replace missing months in reference period according to JRF data
who_dt_totals <- who_dt[!is.na(reference), .(incl = sum(reference), n_incl = .N), by = c("ihme_loc_id", "antigen")]
who_dt_dotals <- who_dt_totals[ihme_loc_id %in% who_dt[is.na(reference) & !is.na(pandemic), ihme_loc_id]]

jrf[, antigen := tolower(Antigen)]
jrf[, Antigen := NULL]

who_dt_totals <- merge(who_dt_totals[n_incl < 12], jrf, by = c("ihme_loc_id", "antigen"), all.x = T)
who_dt_totals[, monthly_fill := (annual_doses - incl)/(12-n_incl)]


who_dt[, ind:= "Raw Value"]
who_dt <- merge(who_dt, who_dt_totals[, .(antigen, ihme_loc_id, monthly_fill)], all.x = T, by = c("antigen", "ihme_loc_id"))
who_dt[is.na(reference), ind := "JRF Imputed"]
who_dt[is.na(reference), reference := monthly_fill]

who_dt <- who_dt[!is.na(reference) & !is.na(pandemic)]

# number plug
unique(who_dt[, .(month_start, location_id)])[, .N, by = "month_start"][order(N)]
who_dt[, sum(reference)+sum(pandemic)]
(who_dt[, .(location_id, month_start)] %>% unique() %>% nrow())*2
who_dt[!(month_start %in% c("jan", "feb"))]

who_dt[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(who_dt)]

who_dt[, age_month_start := 0]
who_dt[, age_month_end := 11]
who_dt[, source_id := 8]
who_dt[, notes := "ratio of ratios, (month 2020/month 2019)/(jan-feb 2020/jan-feb 2019), vaccines administered"]
who_dt[, source := "WHO admin"]

pdf(file.path(out_dir, "plot_who_ratios.pdf"), width = 11, height = 8.5)

for(loc_id in unique(who_dt$location_id)){

  loc_name <- who_dt[location_id == loc_id, unique(location_name)]

  gg <- ggplot(who_dt[location_id == loc_id], aes(x = month))+
    geom_line(aes(y = reference, color = "2019", linetype = ind))+
    geom_line(aes(y = pandemic, color = "2020"))+
    facet_wrap(~antigen)+
    labs(title = loc_name)+
    theme_bw()

  print(gg)


}

who_dt[, ind:= NULL]
who_dt[, monthly_fill := NULL]

dev.off()


# # # Outlier WHO
who_dt <- who_dt[!((location_id == 180 & antigen == "mcv1") | # Kenya MCV1
(location_id == 204) | # Chad
(location_id == 194 & antigen == "mcv1") | # Lesotho MCV
(location_id == 195 & antigen == "mcv1") | # Namibia MCV
(location_id == 216 & antigen == "mcv1") | # Senegal MCV
(location_id == 187 & antigen == "mcv1") | # Somalia MCV
(location_id == 189 & antigen == "mcv1") | # Tanzania MCV
(location_id == 435) | #South Sudan
(location_id == 131) | # Nicaragua
(location_id == 133) | # Venezuela
(location_id == 169 & antigen == "mcv1") | # CAR MCV
(location_id == 17) | # Sri Lanka
(location_id == 18) | # Thailand
(location_id == 218) | # Togo
(location_id == 201))] # Burkina Faso

who_dt <- who_dt[!((location_id == 164) | # Nepal (only use collaborator)
(location_id == 214) | # Nigeria (only use dashboard)
(location_id == 163))] # India (only use dashboard)

# Check for biases in last month of data
who_dt[, max_month := max(month), by = c("location_id", "antigen")]

tmp <- merge(who_dt[month == max_month, .(location_id, location_name, antigen, ratio = pandemic/reference, month)],
             who_dt[month == (max_month-1), .(location_id, antigen, ratio = pandemic/reference)],
             by = c("location_id", "antigen"))

tmp[, change := (ratio.x - ratio.y)/ratio.y]
who_dt <- merge(who_dt, tmp[change<(-0.20) & month>3, .(location_id, antigen, month, flag = 1)],
                all.x = T, by = c("location_id", "antigen", "month"))

tmp[change< (-0.20) & month>3]
tmp[order(change)] %>% head(20)

# # Drop all last month: SENSITIVITY ANALYSIS
# who_dt[month > 3 & month == max_month, flag := 1]

who_dt <- who_dt[is.na(flag)]
who_dt[, c("max_month", "flag") := NULL]

# Nigeria ----------------------------------------------------------------

nga <- dt[source_id ==1 & antigen %in% c("mcv1", "penta3")]

nga[, month_start := substr(month_start, 1, 3)]

nga[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(nga)]

nga[, notes:= "ratio of ratios, (month 2020/month 2019)/(jan-feb 2020/jan-feb 2019), vaccines delivered, administrative"]
nga[, source:= "Nigeria Delivery"]

nga <- nga[, names(who_dt), with = F]

# Nepal -------------------------------------------------------------------

npl <- dt[source_id == 29]

npl[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(npl)]

npl[, notes:= "ratio of ratios, (month 2020/month 2019)/(jan-feb 2020/jan-feb 2019), doses administered"]
npl[, source:= "Nepal admin"]

npl <- npl[, c(names(who_dt)), with = F]


# india -------------------------------------------------------------------

ind <- dt[source_id == 30 & ihme_loc_id == "IND"]

ind[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(ind)]

ind[, notes:= "ratio of ratios, (month 2020/month 2019)/(jan-feb 2020/jan-feb 2019), doses administered"]
ind[, source:= "India HMIS"]

ind <- ind[, c(names(who_dt)), with = F]

ind[is.na(location_id), location_id := gsub("IND_", "", ihme_loc_id)]

# Scotland ----------------------------------------------------------------

sld <- dt[source_id == 25]

sld[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(sld)]
sld[month_start == "july", month:= 7]

sld[, notes:= "proportions comparing monthly rates to pre-pandemic"]
sld[, source:= "PH scotland"]

sld <- sld[, c(names(who_dt)), with = F]


sld[age_month_end == "20w", age_month_end := 5]
sld[, age_month_end := as.integer(age_month_end)]

# eurosurveillance ----------------------------------------------------------------

eng2 <- dt[source_id == 31]

eng2[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(eng2)]

eng2[, notes:= "ratio of ratios, (week 2020/week 2019)/(jan-feb 2020/jan-feb 2019), doses administered, electronic record"]
eng2[, source:= "PH England"]

eng2<- eng2[, c(names(who_dt), "week"), with = F]



# Calculate cumulative ratios for sources with ratio of ratio ---------------------------------------------
ror_dt <- rbindlist(list(who_dt, nga, npl, ind, sld, eng2), fill = T) 
ror_dt[, week := as.integer(week)]

# calculate end date for each source
ror_dt[is.na(week), end_date := timeLastDayInMonth(paste0("2020-", month, "-1")) %>% as.Date()]
ror_dt[!is.na(week), end_date := as.Date("2019-12-28") + 7*week]

pre <- ror_dt[month_start %in% c("jan", "feb"), .(reference = sum(reference), pandemic = sum(pandemic)),
              by = c("location_id", "antigen", "source_id")]

pre[, pre := pandemic/reference]

post <- data.table()

for(e_d in ror_dt[end_date > as.Date("2020-02-29"), unique(end_date)]){
  
  post_tmp <- ror_dt[end_date > as.Date("2020-02-29") & end_date <= e_d, 
                     .(reference = sum(reference), pandemic = sum(pandemic), end_date = as.Date(e_d, origin = "1970-01-01")),
                     by = c("location_id", "antigen", "source_id")]
  
  post <- rbind(post, post_tmp)
}

post[, post := pandemic/reference]

cum_ror_dt <- merge(pre, post, by = c("location_id", "antigen", "source_id"), all = T)
cum_ror_dt[, ratio:= post/pre]
cum_ror_dt[, ratio_se_log := sqrt(1/reference.x + 1/reference.y + 1/pandemic.x + 1/pandemic.y)]
cum_ror_dt[, c("reference.x", "pandemic.x", "reference.y", "pandemic.y") := NULL]

all_ror_dt <- merge(ror_dt[end_date > as.Date("2020-02-29"), 
                       .(location_id, location_name, ihme_loc_id, antigen, end_date, 
                         age_month_start, age_month_end, source_id, source, notes)],
                cum_ror_dt, by = c("location_id", "antigen", "source_id", "end_date"),
                all.x = T)

# MMWR US ----------------------------------------------------------------

usa <- dt[source_id == 4 & age_month_start == 0 & age_month_end == 23] # NOT USING

usa <- usa[, .(location_id, ihme_loc_id, location_name, antigen, source_id, age_month_start, age_month_end, 
                    month_start, week, post = cumsum(pandemic), pre = cumsum(reference),
                    notes = "weekly doses administered 2020/weekly doses administered jan-feb 2020, electronic medical record")]

usa[, ratio:= post/pre]

usa[month_start == "mar" & week == 1, end_date := as.Date("2020-03-07")]
usa[month_start == "mar" & week == 2, end_date := as.Date("2020-03-14")]
usa[month_start == "mar" & week == 3, end_date := as.Date("2020-03-21")]
usa[month_start == "mar" & week == 4, end_date := as.Date("2020-03-28")]
usa[month_start == "mar" & week == 5, end_date := as.Date("2020-04-04")]
usa[month_start == "apr" & week == 1, end_date := as.Date("2020-04-11")]
usa[month_start == "apr" & week == 2, end_date := as.Date("2020-04-18")]
usa[, month_start := NULL]
usa[, week := NULL]
usa[, source := "USA EMR weekly doses"]



# Virgina -----------------------------------------------------------------

vir <- dt[source_id == 24 & antigen %in% c("dtp", "mmr"), .(age_month_start, age_month_end, location_id, ihme_loc_id, location_name,
                                                            source_id,  month_start, week, source = "VA Peds", ratio = cummean(ratio),
                                                            notes = "proportions comparing weekly rates to pre-pandemic"),
          by = "antigen"]


vir[month_start == "mar" & week == 1, end_date := as.Date("2020-03-07")]
vir[month_start == "mar" & week == 2, end_date := as.Date("2020-03-14")]
vir[month_start == "mar" & week == 3, end_date := as.Date("2020-03-21")]
vir[month_start == "mar" & week == 4, end_date := as.Date("2020-03-28")]
vir[month_start == "apr" & week == 1, end_date := as.Date("2020-04-04")]
vir[month_start == "apr" & week == 2, end_date := as.Date("2020-04-11")]

vir[, week := NULL]
vir[, month_start := NULL]


# scientific american ----------------------------------------------------------------

usa_states <- dt[source_id == 5 & location_name != "Washington, D.C."]

usa_states <- usa_states[!(age_month_end %in% c("all", "18yrs"))]

usa_states[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(usa_states)]
usa_states[, end_date := timeLastDayInMonth(paste0("2020-", month, "-1")) %>% as.Date()]
usa_states <- usa_states[order(end_date)]

#extraction error connecticut
usa_states[reference == 1, reference := 100]

pre <- usa_states[month_start %in% c("jan", "feb"), .(reference = mean(reference), pandemic = mean(pandemic)),
           by = c("age_month_start", "age_month_end", "location_id", "ihme_loc_id", "location_name", "antigen", "source_id")]

pre[, pre := pandemic/reference]
pre[, c("reference", "pandemic") := NULL]

post <- usa_states[!(month_start %in% c("jan", "feb")), .(reference = cumsum(reference), pandemic = cumsum(pandemic), end_date),
            by = c("age_month_start", "age_month_end", "location_id", "ihme_loc_id", "location_name", "antigen", "source_id")]

post[, post := pandemic/reference]
post[, c("reference", "pandemic") := NULL]

usa_states <- merge(pre, post, all = T)
usa_states[, ratio:= post/pre]

usa_states[, notes:= "ratio of ratio of ratios, ((month 2020/month 2018)/(month 2019/month 2018))/((jan-feb 2020/jan-feb 2018)/(jan-feb 2019/jan-feb 2018)), state reported doses administered or ordered"]

usa_states[, source := "USA state doses"]

# DC is special, no jan-feb period
dc <- dt[source_id == 5 & location_name == "Washington, D.C."]
dc[location_name == "Washington, D.C.", c("location_id", "ihme_loc_id") := .(531, "USA_531")]

dc[, month:= which(month_start == tolower(month.abb)), by = 1:nrow(dc)]
dc[, end_date := timeLastDayInMonth(paste0("2020-", month, "-1")) %>% as.Date()]
dc <- dc[order(end_date)]

dc <- dc[, .(location_id, ihme_loc_id, location_name, antigen, source_id, age_month_start, age_month_end, 
             post = cumsum(pandemic), pre = cumsum(reference), end_date, source = "USA state doses",
             notes = "ratio of ratios, (month 2020/month 2018)/(month 2019/month 2018), state reported doses administered or ordered")]

dc[, ratio := post/pre]


# ESP ---------------------------------------------------------------------
dt[source_id == 28 & is.na(ratio), ratio:= pandemic/reference]

esp <- dt[source_id == 28, .(age_month_start, age_month_end,
                             location_name, source_id,  month_start, source = "Spain AEV", ratio = cummean(ratio),
                             notes = "proportions comparing monthly rates to same month of 2019"),
          by = c("antigen", "ihme_loc_id")]

esp[, location_id := tstrsplit(ihme_loc_id, "_", keep = 2)]

esp[month_start == "mar", end_date := as.Date("2020-03-31")]
esp[month_start == "apr", end_date := as.Date("2020-04-30")]

esp[, month_start := NULL]

# AUS ---------------------------------------------------------------------

aus <- dt[source_id == 26, 
          .(age_month_start, age_month_end, location_id, ihme_loc_id, location_name,
            source_id, pre = cumsum(reference), post = cumsum(pandemic), month_start, source = "NCIRS AUS",
            notes = "coverage jan-feb 2020 compared to pandemic"), by = "antigen"]

aus[, ratio := post/pre]

aus[, month:= which(substr(month_start, 1, 3) == tolower(month.abb)), by = 1:nrow(aus)]
aus[, end_date := timeLastDayInMonth(paste0("2020-", month, "-1")) %>% as.Date()]

aus[, c("month", "month_start") := NULL]


# Combine -----------------------------------------------------------------

out <- rbindlist(list(all_ror_dt, vir, usa_states, aus, esp, dc), use.names = T, fill = T)
out <- out[order(source_id, ihme_loc_id, end_date)]

out[, age_month := as.integer(age_month_end)]
out[age_month_end == "2yrs", age_month:= 24]

out[, month:= which(end_date == end_week_last_days), by = 1:nrow(out)]
out[is.na(month), month := month(end_date)]

# Only keep monthly data
out <- out[end_date %in% c(last_days, end_week_last_days)]


# replace missing standard error with median by month
median_error <- out[month<=9 & !is.na(ratio_se_log), .(median_error = median(ratio_se_log),
                                            error_lo = quantile(ratio_se_log, .2),
                                            error_hi = quantile(ratio_se_log, .8),
                                            .N), by = "month"]

median_error <- rbind(median_error, 
                          data.table(month = 10:12, 
                                     median_error = median_error[month == 9, median_error],
                                     error_lo = median_error[month == 9, error_lo],
                                     error_hi = median_error[month == 9, error_hi],
                                     N = median_error[month == 9, N]))

ggplot(median_error, aes(x = month, y = median_error))+geom_point()+geom_point(aes(y = error_lo))+geom_point(aes(y = error_hi))
ggplot(median_error, aes(x = month, y = median_error^(-2)))+geom_point()+geom_point(aes(y = error_lo^(-2)))+geom_point(aes(y = error_hi^(-2)))

out <- merge(out, median_error, by = "month", all.x = T)

# replace really small standard error with low percentile by month
out[ratio_se_log < error_lo, ratio_se_log := error_lo]

# replace really large standard error with hi percentile by month
out[ratio_se_log > error_hi, ratio_se_log := error_hi]

# replace missing standard error with median by month
out[is.na(ratio_se_log) & source_id %in% c(26,28), ratio_se_log:= median_error]

# replace missing standard error with high percentile by month for news sources
out[is.na(ratio_se_log), ratio_se_log := error_hi]


out[, c("median_error", "error_hi", "error_lo", "N") := NULL]

out[, weight := 1/ratio_se_log^2]
out[, lower := exp(log(ratio) - 1.96*ratio_se_log)]
out[, upper := exp(log(ratio) + 1.96*ratio_se_log)]
out[, log_ratio := log(ratio)]

write.csv(out, file.path(out_dir, "vaccine_interruption.csv"), row.names = F)

pdf(file.path(out_dir,"vaccine_interruption_scatters.pdf"))

ggplot(out, aes(x=age_month, y=ratio, color = source, size = weight))+geom_point(alpha= 0.5)+scale_y_log10()+scale_x_log10()+theme_bw()
ggplot(out, aes(x=antigen, y=ratio, color = source, size = weight))+geom_point(alpha = 0.5)+scale_y_log10()+theme_bw()
ggplot(out, aes(x = end_date, y = ratio, yend = ratio, color = source), alpha = weight)+geom_point()+scale_y_log10()+theme_bw()

dev.off()

# Set antigens -------------------------------------------------------

out[, c("dtp", "mcv") := 1]
out[antigen %in% c("dtp3", "dpt3", "penta3", "hexavalent", "dtp"), mcv := 0]
out[antigen %in% c("mcv1", "mcv2", "mcv", "mmr", "mmr1"), dtp := 0]


# Impute other antigens -------------------------------------------------

antigen_ratio_dt <- merge(out[dtp == 1 & mcv == 0], 
                          out[mcv == 1 & dtp == 0], 
                          by = c("end_date", "location_id", "source_id"))
antigen_ratio_dt[, weight := 1/(ratio_se_log.x^2 + ratio_se_log.y^2)]
antigen_ratio_dt[, ratio_mcv_dtp := ratio.y/ratio.x]

mcv_dtp_ratio <- antigen_ratio_dt[, weighted.mean(ratio_mcv_dtp, weight)]
# use loc-specific ratio if possible
antigen_ratio_dt_loc <- antigen_ratio_dt[, .(mcv_dtp_ratio = weighted.mean(ratio_mcv_dtp, weight), .N), by = "location_id"]

print(mcv_dtp_ratio)

out_tmp <- data.table()

for(i in 1:nrow(out[dtp==1])){
  
  row <- out[dtp==1][i]
  
  row_mcv <- out[location_id == row$location_id & end_date == row$end_date & source == row$source & mcv == 1]
  
  if(row$location_id %in% antigen_ratio_dt_loc$location_id){
    this_mcv_dtp_ratio <- antigen_ratio_dt_loc[location_id == row$location_id, mcv_dtp_ratio]
  }else{
    this_mcv_dtp_ratio <- mcv_dtp_ratio
  }
  
  if(nrow(row_mcv) ==0){
    tmp_out <- copy(row)
    tmp_out[, ratio := ratio*this_mcv_dtp_ratio]
    tmp_out[, log_ratio := log(ratio)]
    tmp_out[, lower := exp(log_ratio - ratio_se_log*1.96)]
    tmp_out[, upper := exp(log_ratio + ratio_se_log*1.96)]
    tmp_out[, antigen := "mcv imputed"]
    tmp_out[, dtp := 0]
    tmp_out[, mcv := 1]
    
    out_tmp <- rbind(out_tmp, tmp_out)
  }
  
}

for(i in 1:nrow(out[mcv==1])){
  
  row <- out[mcv==1][i]
  
  row_mcv <- out[location_id == row$location_id & end_date == row$end_date & source == row$source & dtp == 1]
  
  if(row$location_id %in% antigen_ratio_dt_loc$location_id){
    this_mcv_dtp_ratio <- antigen_ratio_dt_loc[location_id == row$location_id, mcv_dtp_ratio]
  }else{
    this_mcv_dtp_ratio <- mcv_dtp_ratio
  }
  
  if(nrow(row_mcv) ==0){
    tmp_out <- copy(row)
    tmp_out[, ratio := ratio/this_mcv_dtp_ratio]
    tmp_out[, log_ratio := log(ratio)]
    tmp_out[, lower := exp(log_ratio - ratio_se_log*1.96)]
    tmp_out[, upper := exp(log_ratio + ratio_se_log*1.96)]
    tmp_out[, antigen := "dtp imputed"]
    tmp_out[, dtp := 1]
    tmp_out[, mcv := 0]
    
    out_tmp <- rbind(out_tmp, tmp_out)
  }
  
}

out <- rbind(out, out_tmp)

write.csv(out, file.path(out_dir, "data_mrbrt_source.csv"), row.names = F)

# number plug
out_tmp[, .N, by = "antigen"]
antigen_ratio_dt_loc[location_id %in% out_tmp$location_id]
print(mcv_dtp_ratio)
nrow(antigen_ratio_dt)
unique(antigen_ratio_dt$location_id)

hierarchy[location_id %in% intersect(out_tmp$location_id, antigen_ratio_dt_loc$location_id), 
          .(location_id, location_name)][order(location_name)]

hierarchy[location_id %in% setdiff(out_tmp$location_id, antigen_ratio_dt_loc$location_id), 
          .(location_id, location_name)][order(location_name)]
