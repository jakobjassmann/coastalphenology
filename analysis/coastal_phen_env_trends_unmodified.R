### Coastal Phenology Trend Analysis in Environmental Predictors
### Jakob Assmann j.assmann@ed.ac.uk 16 February 2018

# Dependencies
library(dplyr)
library(MCMCglmm)
library(ggplot2)
library(parallel)
library(gridExtra)

# Housekeeping
script_path <- "scripts/users/jassmann/phenology/coastal_phenology/"
load(file =paste0(script_path, "coastal_phen.Rda"))

# Site names
site_names <- unique(as.character(coastal_phen$site_name))

# Set colours
colour_theme_sites <- data.frame(site_name = c(site_names),
                                 name_pretty =  c("Alexandra Fiord", 
                                                  "Utqiaġvik",
                                                  "Qikiqtaruk",
                                                  "Zackenberg"),
                                 colour = c("#324D5CFF", 
                                            "#46B29DFF", 
                                            "#C2A33EFF",
                                            "#E37B40FF")) # 5th colour #F53855FF

colour_theme_sites$start_year <- c(1989,1993,1999,1994)
colour_theme_sites$end_year <- c(2016,2016,2016,2011)

# Helper function that returns ylab only for Alexfiord (first site)
ylab_filter <- function(label, site_to_plot){
  if(site_to_plot == "ALEXFIORD") {
    return(label)
  } else {
    return("")
  }
}

# Prepare data frame for outputs
trends <- data.frame(site_name = sort(rep(site_names, 3)),
                     predictor = rep(c("Snowmelt", 
                                            "Temperature", 
                                            "Sea-Ice")),
                     intercept = NA,
                     intercept_l95 = NA,
                     intercept_u95 = NA,
                     int_eff_smpl = NA,
                     int_pMCMC = NA,
                     slope = NA,
                     slope_l95 = NA,
                     slope_u95 = NA,
                     slo_eff_smpl = NA,
                     slo_pMCMC = NA)

# Prepare data frame for predictions
years_on_record <- coastal_phen %>%
  select(site_name, year) %>%
  distinct() %>% 
  arrange(site_name,year)
predicts <- data.frame(site_name = rep(years_on_record$site_name, 3),
                       year = as.numeric(rep(years_on_record$year, 3)),
                       predictor = c(rep("Snowmelt", nrow(years_on_record)),
                                     rep("Temperature", nrow(years_on_record)),
                                     rep("Sea-Ice", nrow(years_on_record))),
                       pred = NA,
                       l95 = NA,
                       u95 = NA)

# Define helper function to extract outputs
extr_outputs <- function(site_to_run, predictor){
  if(predictor == "Temperature")
  {model <- get(paste0(site_to_run, "_temp_model"))}
  if(predictor == "Snowmelt")
  {model <- get(paste0(site_to_run, "_snow_model"))}
  if(predictor == "Sea-Ice")
  {model <- get(paste0(site_to_run, "_ice_model"))}
  model_sum <- summary(model)
  trends[trends$site_name == site_to_run &
           trends$predictor == predictor,3:7] <<- model_sum$solutions[1,]
  trends[trends$site_name == site_to_run &
           trends$predictor == predictor,8:12] <<- model_sum$solutions[2,]
}

# Define helper function to calculate predictions
pred_trends <- function(site_to_run, predictor){
  if(predictor == "Temperature"){
    model <- get(paste0(site_to_run, "_temp_model"))
    }
  if(predictor == "Snowmelt"){
    model <- get(paste0(site_to_run, "_snow_model"))
    }
  if(predictor == "Sea-Ice"){
    model <- get(paste0(site_to_run, "_ice_model"))}
  preds <- predict.MCMCglmm(model, interval = "confidence", type = "response")
  preds <- unique(data.frame(preds, year = as.numeric(model$X[,2])))
  preds <- preds[,c(4,1,2,3)]
  if(!identical(predicts[predicts$site_name == site_to_run &
              predicts$predictor == predictor,]$year, preds$year)){
    warning("Years do not match!")
    }
  predicts[predicts$site_name == site_to_run &
             predicts$predictor == predictor, 4:6] <<- preds[,2:4]
}

### site_phen_phase Temperature - Trend & Graph ----

# Definte function to run model for a site
run_temp_model <- function(site_to_run) {
  cphen_sub <- coastal_phen %>%
    filter(site_name == site_to_run) %>%
    select(year, site_phase_temp_unmodified) %>%
    distinct()
  cphen_sub <- as.data.frame(cphen_sub[order(cphen_sub$year),])
  # Parameter extended priors
  model <- MCMCglmm(site_phase_temp_unmodified ~ year,
                    data = cphen_sub,
                    family = "gaussian",
                    nitt = 5000000)
  return(model)  
}
temp_model_list <- lapply(site_names, run_temp_model)
names(temp_model_list) <- make.names(paste0(site_names, "_temp_model"))
list2env(temp_model_list, envir = .GlobalEnv)       
rm(temp_model_list)

# Check convergence
plot(ALEXFIORD_temp_model)
plot(BARROW_temp_model)
plot(QHI_temp_model)
plot(ZACKENBERG_temp_model)

# extract outputs
mapply(extr_outputs, site_names, rep("Temperature",4))

# predictions
mapply(pred_trends, site_names, rep("Temperature", 4))

### Define plotting funciton for a site
plot_site_phase_temp <- function(site_to_plot){
  site_colour <- 
    colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  name_pretty <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$name_pretty
  
  site_phen <- coastal_phen %>%
    filter(site_name == site_to_plot) %>%
    select(year, site_phase_temp_unmodified) %>%
    distinct()
  site_phen <- as.data.frame(site_phen[order(site_phen$year),])
  
  preds <- predicts %>% 
    filter(site_name == site_to_plot, predictor == "Temperature")
    
  ggplot(data = site_phen, aes(x = year, y = site_phase_temp_unmodified)) +
    geom_point(size = 4, colour = site_colour) +
    geom_line(data = preds, aes(x = year, y = pred),
              colour = site_colour,
              inherit.aes = F,
              size = 1) +
    geom_ribbon(data = preds, aes(x = year, ymin = l95, ymax = u95),
                fill = site_colour,
                alpha = 0.5, 
                inherit.aes = F) +
    scale_x_continuous(limits= c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year, 
                                 colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year), 
                       breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(-5.2, 10), breaks = seq(-5, 10, 5)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Spring Temperature (C)", site_to_plot)) +
    annotate("text", x = (colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year + 
                            colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year) / 2,
             y = 9.7, label = name_pretty, 
             colour = site_colour, size = 7, hjust = 0.5) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 20),
          legend.position = "none")
}

# Make plots for all sites
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_site_phase_temp_plot"))),
                plot_site_phase_temp), 
         envir = .GlobalEnv)

# Create panel and save to file
panel_layout <- rbind(c(1,2,3,4))

plot_grob <- grid.arrange(ALEXFIORD_site_phase_temp_plot, 
                          BARROW_site_phase_temp_plot, 
                          QHI_site_phase_temp_plot, 
                          ZACKENBERG_site_phase_temp_plot,
                          layout_matrix = panel_layout)

ggsave(filename = paste0(script_path, "/coastal_site_phase_temp_plot_MCMC.png"), 
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)

### 

### Snowmelt - Trend & Graph ----
# Definte function to run model for a site
run_snow_model <- function(site_to_run) {
  cphen_sub <- coastal_phen %>%
    filter(site_name == site_to_run) %>%
    group_by(year) %>%
    summarise(snowmelt = mean(snowmelt, na.rm = T))
  cphen_sub <- as.data.frame(cphen_sub[order(cphen_sub$year),])
  model <- MCMCglmm(snowmelt ~ year, data = cphen_sub,
                    family = "gaussian",
                    nitt = 5000000)
  return(model)  
}
snow_model_list <- lapply(site_names, run_snow_model)
names(snow_model_list) <- make.names(paste0(site_names, "_snow_model"))
list2env(snow_model_list, envir = .GlobalEnv)       
rm(snow_model_list)

# Check convergence
plot(ALEXFIORD_snow_model)
plot(BARROW_snow_model)
plot(QHI_snow_model)
plot(ZACKENBERG_snow_model)

# extract outputs
mapply(extr_outputs, site_names, rep("Snowmelt",4))

# predictions
mapply(pred_trends, site_names, rep("Snowmelt", 4))

# Define snowmelt plot function
plot_snowmelt <- function(site_to_plot){
  site_colour <- 
    colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  name_pretty <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$name_pretty
  
  site_phen <- coastal_phen %>% 
    filter(site_name == site_to_plot) %>% 
    group_by(year) %>% 
    summarise(snowmelt = round(mean(snowmelt, na.rm = T)))

  preds <- predicts %>% 
    filter(site_name == site_to_plot, predictor == "Snowmelt")
  
  ggplot(data = site_phen, aes(x = year, y = snowmelt)) +
    geom_point(size = 4, colour = site_colour) +    
    geom_line(data = preds, aes(x = year, y = pred),
              colour = site_colour,
              inherit.aes = F,
              size = 1) +
    geom_ribbon(data = preds, aes(x = year, ymin = l95, ymax = u95),
                fill = site_colour,
                alpha = 0.5, 
                inherit.aes = F) +
    scale_x_continuous(limits= c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year, 
                                 colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year),
                       breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(120, 200), breaks = seq(100, 240, 20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Snowmelt Date (DoY)", site_to_plot)) +
    annotate("text", x = (colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year + 
                            colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year) / 2,
             y = 200, label = name_pretty, 
             colour = site_colour, size = 7, hjust = 0.5) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 20),
          legend.position = "none")
}
# apply to all sites
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_snowmelt_plot"))),
                plot_snowmelt), 
         envir = .GlobalEnv)
# plot
panel_layout <- rbind(c(1,2,3,4))

plot_grob <- grid.arrange(ALEXFIORD_snowmelt_plot, 
                          BARROW_snowmelt_plot, 
                          QHI_snowmelt_plot, 
                          ZACKENBERG_snowmelt_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_snowmelt_plot_MCMC.png"),
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)

### Sea-Ice - Trend & Graph ----

# Definte function to run model for a site
run_ice_model <- function(site_to_run) {
  cphen_sub <- coastal_phen %>%
    filter(site_name == site_to_run) %>%
    select(year, onset_ice_melt) %>%
    distinct()
  cphen_sub <- as.data.frame(cphen_sub[order(cphen_sub$year),])
  model <- MCMCglmm(onset_ice_melt ~ year, data = cphen_sub,
                    family = "gaussian",
                    nitt = 5000000)
  return(model)  
}
ice_model_list <- lapply(site_names, run_ice_model)
names(ice_model_list) <- make.names(paste0(site_names, "_ice_model"))
list2env(ice_model_list, envir = .GlobalEnv)       
rm(ice_model_list)

# Check convergence
plot(ALEXFIORD_ice_model)
plot(BARROW_ice_model)
plot(QHI_ice_model)
plot(ZACKENBERG_ice_model)

# extract outputs
mapply(extr_outputs, site_names, rep("Sea-Ice",4))

# predictions
mapply(pred_trends, site_names, rep("Sea-Ice", 4))

# Define sea ice  plot function
plot_sea_ice <- function(site_to_plot){
  site_colour <- 
    colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  name_pretty <- 
    colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$name_pretty
  
  site_phen <- coastal_phen %>%
    filter(site_name == site_to_plot) %>%
    select(year, onset_ice_melt) %>%
    distinct()
  
  preds <- predicts %>% 
    filter(site_name == site_to_plot, predictor == "Sea-Ice")
  
  ggplot(data = site_phen, aes(x = year, y = onset_ice_melt)) +
    geom_point(size = 4, colour = site_colour) +
    geom_line(data = preds, aes(x = year, y = pred),
              colour = site_colour,
              inherit.aes = F,
              size = 1) +
    geom_ribbon(data = preds, aes(x = year, ymin = l95, ymax = u95),
                fill = site_colour,
                alpha = 0.5, 
                inherit.aes = F) +
    scale_x_continuous(limits= c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year, 
                                         colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year),
                       breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(90, 250), breaks = seq(100,240,20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Spring Drop Sea-Ice (DoY)", site_to_plot)) +
    annotate("text", x = (colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year + 
                            colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year) / 2,
             y = 245, label = name_pretty,
             colour = site_colour, size = 7, hjust = 0.5) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 20),
          legend.position = "none")
}
# apply to all sites
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_sea_ice_plot"))),
                plot_sea_ice), 
         envir = .GlobalEnv)
# plot
panel_layout <- rbind(c(1,2,3,4))

plot_grob <- grid.arrange(ALEXFIORD_sea_ice_plot, 
                          BARROW_sea_ice_plot, 
                          QHI_sea_ice_plot, 
                          ZACKENBERG_sea_ice_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_sea_ice_plot_MCMC.png"),
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)


### Export trends to file ----

# Change per decade 
trends <- trends %>% mutate(adv_per_dec = slope * 10)

write.csv(trends, file = paste0(script_path, "coastal_phen_env_trends.csv"), 
          row.names = F)

coastal_phen %>% group_by(site_name) %>% summarise(mean_interval = mean((doy - prior_visit), na.rm = T))

### Export models to files ----
model_names <- ls()[grep("_model", ls())]
model_names <- model_names[-grep("run", model_names)]
save(list = model_names, file = paste0(script_path, "models/env_models.Rda"))

# load if previously run
load(file = paste0(script_path, "models/env_models.Rda"))
