#-------------------Header------------------------------------------------
# Author: NAME
# Date: 2/8/21
# Purpose: Make paper plots
#          
# source("FILEPATH/paper_plots.R", echo=T)
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

packages <- c("data.table","magrittr","ggplot2", "lubridate", "gridExtra", "grid", "stringr", "scales", "cowplot")

for(p in packages){
  if(p %in% rownames(installed.packages())==FALSE){
    install.packages(p)
  }
  library(p, character.only = T)
}

# Directories -------------------------------------------------------------

model_name <- "VERSION_NAME"

work_dir <- file.path(j_root, "FILEPATH", paste0("cumulative", model_name))
output_dir <- file.path(j_root, "FILEPATH", model_name)
dir.create(output_dir)


source(file.path(central_lib, "FILEPATH/get_location_metadata.R"))
hierarchy <- get_location_metadata(35, gbd_round_id = 6)


model_names <- c("vaccine_dtp", "vaccine_mcv")

color_scheme <- c("DTP3" = "#a8216b",
                  "MCV1" = "#57c1e4")

color_scheme_reg <- c("Central Europe, Eastern Europe, and Central Asia" = "#0A2F51",
                      "High-income" = "#57c1e4",
                      "Latin America and Caribbean" = "#2e9598",
                      "North Africa and Middle East" = "#f7db69",
                      "Southeast Asia, East Asia, and Oceania" = "#f26a44",
                      "South Asia" = "#ec1b4b",
                      "Sub-Saharan Africa" = "#a8216b")




# Disruption line plot ----------------------------------------------------
dt <- fread(file.path(output_dir, "summary_results_agg.csv"))

plot <- dt[period %in% month.abb & level <=1]
plot[, period := factor(period, levels = month.abb)]
plot[, month := as.numeric(period)]
plot[model == "vaccine_dtp", vaccine := "DTP3"]
plot[model == "vaccine_mcv", vaccine := "MCV1"]


plot[, location_name := factor(location_name, levels = c("Global", sort(unique(dt[level == 1, location_name]))))]

relabel <- function(string){
  paste(strwrap(string, width = 40),collapse="\n")
}

loc_labs <- lapply(levels(plot$location_name), relabel) %>% unlist
names(loc_labs) <- levels(plot$location_name)

pdf(file = paste0(output_dir, "/regional_coverage_lineplot.pdf"),
    onefile = TRUE,
    width = 5.5,
    height = 10)

gg <- ggplot(plot, aes(x = month, y = mean_coverage * 100, ymin = lower_coverage *100, ymax = upper_coverage*100,
                       color = vaccine, fill = vaccine))+
  geom_hline(aes(yintercept = mean_coverage_ref*100, linetype = "Expected coverage", color = vaccine))+
  geom_line(aes(linetype = "Estimated coverage"))+
  geom_ribbon(alpha = 0.1, color = NA, show.legend = F)+
  scale_color_manual(values = color_scheme)+
  scale_fill_manual(values = color_scheme)+
  facet_wrap( ~ location_name, labeller = labeller(location_name = loc_labs), ncol = 2)+
  scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), 
                     expand = c(0,0), minor_breaks = c(1:12))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 120), breaks = c(20, 40, 60, 80, 100))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  labs(x = "Month", y = "Vaccine coverage (%)", color = "", linetype = "")+
  guides(fill = FALSE)+theme_bw()+
  theme(legend.position = "bottom", panel.spacing = unit(0.5, "lines"))


print(gg)

dev.off()




# number of missed vaccines plot--------------------------------------------------------------


pdf(file = paste0(output_dir, "/regional_num_vaccines_missed.pdf"),
    onefile = TRUE,
    width = 11,
    height = 6.5)


plot <- dt[period %in% month.abb & level == 1]
plot[, period := factor(period, levels = month.abb)]
plot[, month := as.numeric(period)]
plot[model == "vaccine_dtp", vaccine := "DTP3"]
plot[model == "vaccine_mcv", vaccine := "MCV1"]

plot[, location_name := factor(location_name, levels = names(color_scheme_reg))]

test <- melt.data.table(plot, id.vars = c("location_name", "period", "month", "vaccine"), measure.vars = c("mean_n_missed", "mean_n_missed_ref"))

test[variable == "mean_n_missed_ref", glob_val := sum(value), by = c("month", "vaccine")]

gg_glob <- ggplot(test[variable == "mean_n_missed"], aes(x= month, y = value/1e6, fill = location_name))+
  geom_col()+
  geom_hline(data = test[variable == "mean_n_missed_ref"], aes(yintercept = glob_val/1e6), color = "white")+
  geom_hline(data = test[variable == "mean_n_missed_ref"], aes(yintercept = glob_val/1e6), linetype = "dashed")+
  facet_wrap(~vaccine)+
  scale_fill_manual(values = color_scheme_reg)+
  scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), 
                     expand = c(0,0))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  labs(title = "A)", x = "Month", y = "Number of children missing dose (millions)")+
  theme_bw()+
  theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "lines"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "lines"))


gg_reg <- ggplot(test[variable == "mean_n_missed"], aes(x= month, y = value/1e3, fill = location_name))+
  geom_col()+
  geom_hline(data = test[variable == "mean_n_missed_ref"], aes(yintercept = value/1e3), color = "white")+
  geom_hline(data = test[variable == "mean_n_missed_ref"], aes(yintercept = value/1e3, 
                                                               linetype = "Expected number of children missing doses in the absence of COVID-19"))+
  facet_grid(location_name~vaccine, scales = "free_y")+
  scale_fill_manual(values = color_scheme_reg,
                    labels = wrap_format(30))+
  scale_x_continuous(breaks = c(1:12), labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"), 
                     expand = c(0,0))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_linetype_manual(values = "dashed",
                        labels = wrap_format(30))+
  labs(title = "B)", x = "Month", y = "Number of children missing dose (thousands)", fill = "", linetype = "")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.spacing = unit(0.5, "lines"), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "lines"),
        strip.background.y = element_blank(), strip.text.y = element_blank(),
        legend.text = element_text(size = 8))


p_all <- plot_grid(plotlist = list(gg_glob, gg_reg), ncol = 2, rel_widths = c(0.4, 0.6))

print(p_all)


dev.off()
