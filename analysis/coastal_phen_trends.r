### Coastal Phenology Trend Analysis in Phenology Observations and Environmental
### Predictors by Jakob Assmann 13 March 2018

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

# Site spp phen combos
site_spp_phen_combos <- unique(coastal_phen$site_spp_phen)

# Spp phen combos
spp_phen_combos <- sort(unique(coastal_phen$site_spp_phen))

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

colour_theme_spp_phen <- data.frame(site_name = c(rep(site_names[1], 8), 
                                                  rep(site_names[2], 7),
                                                  rep(site_names[3], 3),
                                                  rep(site_names[4], 6)),
                                    site_spp_phen = spp_phen_combos,
                                    spp_phen = c(unique(sort(as.character(coastal_phen[coastal_phen$site_name == site_names[1],]$spp_phen))),
                                                 unique(sort(as.character(coastal_phen[coastal_phen$site_name == site_names[2],]$spp_phen))),
                                                 unique(sort(as.character(coastal_phen[coastal_phen$site_name == site_names[3],]$spp_phen))),
                                                 unique(sort(as.character(coastal_phen[coastal_phen$site_name == site_names[4],]$spp_phen)))),
                                    colour = c(
                                      colorRampPalette(c("#19262EFF","#84CBF2FF"), 
                                                       interpolate = "linear")(8), # ALEXFIIORD
                                      colorRampPalette(c("#112B26FF","#5BE8CDFF"), 
                                                       interpolate = "linear")(7), # BARROW
                                      colorRampPalette(c("#382F12FF","#C2A33EFF"), 
                                                       interpolate = "linear")(3), # QHI
                                      colorRampPalette(c("#422413FF","#FF8A48FF"), 
                                                       interpolate = "linear")(6) # ZACKENBERG
                                    ), stringsAsFactors = F)

# Helper function that returns ylab only for Alexfiord (first site)
ylab_filter <- function(label, site_to_plot){
  if(site_to_plot == "ALEXFIORD") {
    return(label)
  } else {
    return("")
  }
}

### Quality Control Phenology data----

## Check for NAs and quality in predictors and responses

# Phenology observations, upper and lower bounds
unique(coastal_phen$doy)
unique(coastal_phen$prior_visit)
sum(is.na(coastal_phen$doy))
sum(is.na(coastal_phen$prior_visit))
# remove NAs (211 values)
coastal_phen <- coastal_phen[!is.na(coastal_phen$doy),]
# nice!

# Snowmelt dates
unique(coastal_phen$snowmelt)
sum(is.na(coastal_phen$snowmelt))
# Some NAs in there (153)
coastal_phen <- coastal_phen[!is.na(coastal_phen$snowmelt),]
# Nice!

# Remaining predictors
unique(coastal_phen$phase_temp)
unique(coastal_phen$spring_temp)
unique(coastal_phen$apr_temp)
unique(coastal_phen$may_temp)
unique(coastal_phen$jun_temp)
unique(coastal_phen$onset_ice_melt)
unique(coastal_phen$spring_extent)

sum(is.na(coastal_phen$phase_temp))
sum(is.na(coastal_phen$spring_temp))
sum(is.na(coastal_phen$apr_temp))
sum(is.na(coastal_phen$may_temp))
sum(is.na(coastal_phen$jun_temp))
sum(is.na(coastal_phen$onset_ice_melt))
sum(is.na(coastal_phen$spring_extent))
# Okay. All good.


### Trend analysis fro phenological observations ----

# If previously run load from this file:
# load(paste0(script_path, "models/phen_trend_models.Rda"))

# Define function to run MCMCglmm models
trend_model <- function(site_spp_sphen_select){
  cphen_subset <- coastal_phen %>% 
    filter(site_spp_phen == site_spp_sphen_select) 
    model <- MCMCglmm(cbind(prior_visit, doy) ~ year, 
                      random = ~ year_fac + plot_id,
                      data = as.data.frame(cphen_subset),
                      family = "cengaussian",
                      nitt = 1000000)
    return(model)
}

# Execute over all site_spp_phen_combinations with parallel processing
# models <-  mclapply(site_spp_phen_combos, trend_model, mc.cores = 4)
# names(models) <- make.names(paste0(site_spp_phen_combos, "_model"))
# list2env(models, envir = .GlobalEnv)

# Dont'r run in parallel, as we will not know which model will be associated
# with which output afterwards. Better run linearly.
model_list <- lapply(site_spp_phen_combos, trend_model)
names(model_list) <- make.names(paste0(site_spp_phen_combos, "_model"))
list2env(model_list, envir = .GlobalEnv)                                

### Gather outputs 
## Prep data frame
phen_trends <- data.frame(site_spp_phen = site_spp_phen_combos,
                          intercept = NA,
                          intercept_l95 = NA,
                          intercept_u95 = NA,          
                          intercept_esmpl = NA,
                          intercept_pMCMC = NA,
                          slope = NA,
                          slope_l95 = NA,
                          slope_u95 = NA,
                          slope_esmpl = NA,
                          slope_pMCMC = NA)
## Define funciton to extract slope parameters and estimates 
extr_outputs <- function(site_spp_sphen_select){
  model <- get(paste0(site_spp_sphen_select, "_model"))
  model_sum <- summary(model)
  # retrive row index for output dataframe
  row_no <- which(phen_trends$site_spp_phen == site_spp_sphen_select)
  # use one liners to extract intersept and slope estimators
  phen_trends[row_no,2:6] <<- model_sum$solutions[1,]
  phen_trends[row_no,7:11] <<- model_sum$solutions[2,]
}

## Apply to all site_spp_phen combos
lapply(site_spp_phen_combos, extr_outputs)

### Claculate predictions
# prep data frame
phen_predicts <- data.frame(site_spp_phen = coastal_phen %>% 
                              select(site_spp_phen, year) %>%
                              distinct(site_spp_phen, year) %>%
                              select(site_spp_phen) %>%
                              unlist() %>% 
                              as.character(),
                            year = coastal_phen %>% 
                              select(site_spp_phen, year) %>%
                              distinct(site_spp_phen, year) %>%
                              select(year) %>%
                              unlist() %>% 
                              as.numeric(),
                            pred = NA,
                            pred_l95 = NA,
                            pred_u95 = NA)

# Define function to extract predicitons
extract_pred <- function(site_spp_phen_select){
  model <- get(paste0(site_spp_phen_select, "_model"))
  row_nos <- which(phen_predicts$site_spp_phen == site_spp_phen_select)
  preds <- predict.MCMCglmm(model, interval = "confidence", type = "response")
  preds <- unique(data.frame(preds, year = as.numeric(model$X[,2])))
  preds <- preds[,c(4,1,2,3)]
  phen_predicts[row_nos, 3:5] <<- preds[,2:4]
  paste0("Done with: ", site_spp_phen_select, " ;")
}
lapply(site_spp_phen_combos, extract_pred)

# create site name and spp_phen columns
phen_predicts$site_name <- gsub("_.*", "", phen_predicts$site_spp_phen) 
phen_predicts$spp_phen <- gsub("^[A-Z]*_", "", phen_predicts$site_spp_phen) 

### Plot phenology data with trends ----
# First without ZACK SILACA due to large errors
colour_theme_sites$start_year <- c(1989,1993,1999,1994)
colour_theme_sites$end_year <- c(2016,2016,2016,2011)
colour_theme_sites$early_doy <- c(130,130,130,130)
colour_theme_sites$late_doy <- c(220,220,220,220)

plot_phen <- function(site_to_plot){

  site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  name_pretty <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$name_pretty
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(spp_phen, year) %>% mutate(int_mean = (prior_visit + doy) /2) %>% 
    summarise(phen_mean = round(mean(int_mean, na.rm = T)))
  site_phen$spp_phen <- factor(as.character(site_phen$spp_phen))
  site_preds <- phen_predicts %>%
    filter(site_name == site_to_plot)
  site_preds$spp_phen <- factor(site_preds$spp_phen)
  ribbon_alpha <- 0.5
  site_colours_ribbon <- site_colours
  if(site_to_plot == "ZACKENBERG"){
    site_phen$spp_phen <- factor(site_phen$spp_phen, levels = rev(levels(site_phen$spp_phen)))
    site_preds$spp_phen <- factor(site_preds$spp_phen, levels = rev(levels(site_preds$spp_phen)))
    site_colours <- rev(site_colours)
    # Exclude SILACA due to large error bars
    site_colours_ribbon <- c("#FF8A4800", "#D9753D80", "#B3613280", "#8D4C2880", "#67381D80", "#42241380")
    
    ribbon_alpha <- NA
    
  }
  
  ggplot() +
    geom_ribbon(data = site_preds[,], #-which(site_preds$spp_phen == "SILACA_flowering")
                mapping = aes(x = year,
                              ymin = pred_l95,
                              ymax = pred_u95,
                              fill = spp_phen,
                              group = spp_phen),
                alpha = ribbon_alpha,
                inherit.aes = FALSE) +
    geom_point(data = site_phen, 
               mapping = aes(x = year, 
                             y = phen_mean, 
                             colour = spp_phen, 
                             fill= spp_phen, 
                             group = spp_phen),
               size = 4) +
    geom_line(data = site_preds,
              mapping = aes(x = year,
                            y = pred,
                            colour = spp_phen,
                            group = spp_phen),
              inherit.aes = FALSE,
              size = 1) +
    scale_colour_manual(values = site_colours) +
    scale_fill_manual(values = site_colours_ribbon) +

    scale_x_continuous(limits= c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year, 
                                 colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year), 
                       breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$early_doy, 
                                  colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$late_doy), 
                       breaks = seq(100, 240, 20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Mean Phenology (DoY)", site_to_plot)) +
    annotate("text", x = (colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year + 
               colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year) / 2, 
             y = colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$late_doy, 
             label = name_pretty, 
             colour = site_colour, size = 7, hjust = 0.5) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 15),
          legend.position = "none") 
}
# plot all sites
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_phen_plot"))),
                plot_phen), 
         envir = .GlobalEnv)

# Set panel layout
panel_layout <- rbind(c(1,2,3,4))

# Combine plots and export
plot_grob <- grid.arrange(ALEXFIORD_phen_plot, 
                          BARROW_phen_plot, 
                          QHI_phen_plot, 
                          ZACKENBERG_phen_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_phen_plot_MCMCglmm_intmeans.png"), 
       plot = plot_grob, scale = 1.7, 
       width = 10, 
       height = (10/3), 
       dpi = 600)

# Second with ZACK SILACA due to large errors
plot_phen <- function(site_to_plot){
  
  site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  name_pretty <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$name_pretty
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(spp_phen, year) %>% mutate(int_mean = (prior_visit + doy) /2) %>% 
    summarise(phen_mean = round(mean(int_mean, na.rm = T)))
  site_phen$spp_phen <- factor(as.character(site_phen$spp_phen))
  site_preds <- phen_predicts %>%
    filter(site_name == site_to_plot)
  site_preds$spp_phen <- factor(site_preds$spp_phen)
  ribbon_alpha <- 0.5
  site_colours_ribbon <- site_colours
  if(site_to_plot == "ZACKENBERG"){
    site_phen$spp_phen <- factor(site_phen$spp_phen, levels = rev(levels(site_phen$spp_phen)))
    site_preds$spp_phen <- factor(site_preds$spp_phen, levels = rev(levels(site_preds$spp_phen)))
    site_colours <- rev(site_colours)
    # Exclude SILACA due to large error bars
    site_colours_ribbon <- c("#FF8A4880", "#D9753D80", "#B3613280", "#8D4C2880", "#67381D80", "#42241380")
    
    ribbon_alpha <- NA
  }
  
  ggplot() +
    geom_ribbon(data = site_preds[,], #-which(site_preds$spp_phen == "SILACA_flowering")
                mapping = aes(x = year,
                              ymin = pred_l95,
                              ymax = pred_u95,
                              fill = spp_phen,
                              group = spp_phen),
                alpha = ribbon_alpha,
                inherit.aes = FALSE) +
    geom_point(data = site_phen, 
               mapping = aes(x = year, 
                             y = phen_mean, 
                             colour = spp_phen, 
                             fill= spp_phen, 
                             group = spp_phen),
               size = 4) +
    geom_line(data = site_preds,
              mapping = aes(x = year,
                            y = pred,
                            colour = spp_phen,
                            group = spp_phen),
              inherit.aes = FALSE,
              size = 1) +
    scale_colour_manual(values = site_colours) +
    scale_fill_manual(values = site_colours_ribbon) +
    
    scale_x_continuous(limits= c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year, 
                                 colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year), 
                       breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$early_doy, 
                                  colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$late_doy), 
                       breaks = seq(100, 240, 20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Mean Phenology (DoY)", site_to_plot)) +
    annotate("text", x = (colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$end_year + 
                            colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$start_year) / 2, 
             y = colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$late_doy, 
             label = name_pretty, 
             colour = site_colour, size = 7, hjust = 0.5) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 15),
          legend.position = "none") 
}
# plot all sites
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_phen_plot"))),
                plot_phen), 
         envir = .GlobalEnv)

# Set panel layout
panel_layout <- rbind(c(1,2,3,4))

# Combine plots and export
plot_grob <- grid.arrange(ALEXFIORD_phen_plot, 
                          BARROW_phen_plot, 
                          QHI_phen_plot, 
                          ZACKENBERG_phen_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_phen_plot_MCMCglmm_intmeans_with_SILACA.png"), 
       plot = plot_grob, scale = 1.7, 
       width = 10, 
       height = (10/3), 
       dpi = 600)
### Calculate decadal advance ---- 
phen_trends <- phen_trends %>% mutate(adv_per_dec = slope * 10)

# Min and max advance
phen_trends %>% summarise(min_adv_per_dec = round(min(adv_per_dec),2), 
                          max_adv_per_dec = round(max(adv_per_dec),2))

# Export
write.csv(phen_trends, file = paste0(script_path, "coastal_phen_trends.csv"), row.names = F)


# And safe the model files for the record

model_names <- ls()[grep("_model", ls())][-19] # excluding the model function
save(list = model_names, file = paste0(script_path, "models/phen_trend_models.Rda"))

# Check Convergance
lapply(model_names, function(x) { 
  
  png(file = paste0(script_path, "../quality_control/phen_trends/", x, "_sol.png"))
  par(ask = F)
  plot(get(x)$Sol)
  dev.off()
  
  png(file = paste0(script_path, "../quality_control/phen_trends/", x, "_vcv.png"))
  par(ask = F)
  plot(get(x)$VCV)
  dev.off()
  
  }
  )
# Check summaries
for(x in 1:length(model_names)){
  print(model_names[x])
  print(summary(get(model_names[x])))
  readline("Hit button to see next summary:")
  }
summary(get(model_names[1]))

# 
# ####### Aditional calcs: 1) Begining of season data -5 instead of -10 ----
# 
# rm(list = ls())
# # Housekeeping
# script_path <- "scripts/users/jassmann/phenology/coastal_phenology/"
# load(file =paste0(script_path, "coastal_phen.Rda"))
# 
# # remove NAs (211 values)
# coastal_phen <- coastal_phen[!is.na(coastal_phen$doy),]
# 
# # Site names
# site_names <- unique(as.character(coastal_phen$site_name))
# 
# # Site spp phen combos
# site_spp_phen_combos <- unique(coastal_phen$site_spp_phen)
# 
# # Spp phen combos
# spp_phen_combos <- sort(unique(coastal_phen$site_spp_phen))
# 
# # Set colours
# colour_theme_sites <- data.frame(site_name = c(site_names),
#                                  colour = c("#324D5CFF", 
#                                             "#46B29DFF", 
#                                             "#C2A33EFF",
#                                             "#E37B40FF")) # 5th colour #F53855FF
# colour_theme_spp_phen <- data.frame(site_name = c(rep(site_names[1], 8), 
#                                                   rep(site_names[2], 7),
#                                                   rep(site_names[3], 3),
#                                                   rep(site_names[4], 6)),
#                                     site_spp_phen = spp_phen_combos,
#                                     spp_phen = c(unique(as.character(coastal_phen[coastal_phen$site_name == site_names[1],]$spp_phen)),
#                                                  unique(as.character(coastal_phen[coastal_phen$site_name == site_names[2],]$spp_phen)),
#                                                  unique(as.character(coastal_phen[coastal_phen$site_name == site_names[3],]$spp_phen)),
#                                                  unique(as.character(coastal_phen[coastal_phen$site_name == site_names[4],]$spp_phen))),
#                                     colour = c(
#                                       colorRampPalette(c("#19262EFF","#84CBF2FF"), 
#                                                        interpolate = "linear")(8), # ALEXFIIORD
#                                       colorRampPalette(c("#112B26FF","#5BE8CDFF"), 
#                                                        interpolate = "linear")(7), # BARROW
#                                       colorRampPalette(c("#382F12FF","#C2A33EFF"), 
#                                                        interpolate = "linear")(3), # QHI
#                                       colorRampPalette(c("#422413FF","#FF8A48FF"), 
#                                                        interpolate = "linear")(6) # ZACKENBERG
#                                     ), stringsAsFactors = F)
# 
# # Helper function that returns ylab only for Alexfiord (first site)
# ylab_filter <- function(label, site_to_plot){
#   if(site_to_plot == "ALEXFIORD") {
#     return(label)
#   } else {
#     return("")
#   }
# }
# 
# 
# # adjust data set
# coastal_phen[(coastal_phen$doy - coastal_phen$prior_visit) == 10,
#              ]$prior_visit <- 
#   coastal_phen[(coastal_phen$doy - coastal_phen$prior_visit) == 10,
#                ]$doy - 5
# # Define function to run MCMCglmm models
# trend_model <- function(site_spp_sphen_select){
#   cphen_subset <- coastal_phen %>% 
#     filter(site_spp_phen == site_spp_sphen_select) 
#   model <- MCMCglmm(cbind(prior_visit, doy) ~ year, data = as.data.frame(cphen_subset),
#                     family = "cengaussian",
#                     nitt = 20000)
#   return(model)
# }
# 
# # Execute over all site_spp_phen_combinations
# model_list <- lapply(site_spp_phen_combos, trend_model)
# names(model_list) <- make.names(paste0(site_spp_phen_combos, "_model"))
# list2env(model_list, envir = .GlobalEnv)                                
# 
# ### Gather outputs 
# ## Prep data frame
# phen_trends <- data.frame(site_spp_phen = site_spp_phen_combos,
#                           intercept = NA,
#                           intercept_l95 = NA,
#                           intercept_u95 = NA,          
#                           intercept_esmpl = NA,
#                           intercept_pMCMC = NA,
#                           slope = NA,
#                           slope_l95 = NA,
#                           slope_u95 = NA,
#                           slope_esmpl = NA,
#                           slope_pMCMC = NA)
# ## Define funciton to extract slope parameters and estimates 
# extr_outputs <- function(site_spp_sphen_select){
#   model <- get(paste0(site_spp_sphen_select, "_model"))
#   model_sum <- summary(model)
#   # retrive row index for output dataframe
#   row_no <- which(phen_trends$site_spp_phen == site_spp_sphen_select)
#   # use one liners to extract intersept and slope estimators
#   phen_trends[row_no,2:6] <<- model_sum$solutions[1,]
#   phen_trends[row_no,7:11] <<- model_sum$solutions[2,]
# }
# 
# ## Apply to all site_spp_phen combos
# lapply(site_spp_phen_combos, extr_outputs)
# 
# ### Claculate predictions
# # prep data frame
# phen_predicts <- data.frame(site_spp_phen = coastal_phen %>% 
#                               select(site_spp_phen, year) %>%
#                               distinct(site_spp_phen, year) %>%
#                               select(site_spp_phen) %>%
#                               unlist() %>% 
#                               as.character(),
#                             year = coastal_phen %>% 
#                               select(site_spp_phen, year) %>%
#                               distinct(site_spp_phen, year) %>%
#                               select(year) %>%
#                               unlist() %>% 
#                               as.numeric(),
#                             pred = NA,
#                             pred_l95 = NA,
#                             pred_u95 = NA)
# 
# # Define function to extract predicitons
# extract_pred <- function(site_spp_phen_select){
#   model <- get(paste0(site_spp_phen_select, "_model"))
#   row_nos <- which(phen_predicts$site_spp_phen == site_spp_phen_select)
#   preds <- predict.MCMCglmm(model, interval = "confidence", type = "response")
#   preds <- unique(data.frame(preds, year = as.numeric(model$X[,2])))
#   preds <- preds[,c(4,1,2,3)]
#   phen_predicts[row_nos, 3:5] <<- preds[,2:4]
#   paste0("Done with: ", site_spp_phen_select, " ;")
# }
# lapply(site_spp_phen_combos, extract_pred)
# 
# # create site name and spp_phen columns
# phen_predicts$site_name <- gsub("_.*", "", phen_predicts$site_spp_phen) 
# phen_predicts$spp_phen <- gsub("^[A-Z]*_", "", phen_predicts$site_spp_phen) 
# 
# ### Plot phenology data with trends ----
# plot_phen <- function(site_to_plot){
#   
#   site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
#   site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
#   site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
#     group_by(spp_phen, year) %>% summarise(phen_mean = round(mean(doy, na.rm = T)))
#   
#   site_preds <- phen_predicts %>%
#     filter(site_name == site_to_plot)
#   
#   ggplot(data = site_phen, aes(x = year, 
#                                y = phen_mean, 
#                                colour = spp_phen, 
#                                fill= spp_phen, 
#                                group = spp_phen)) +
#     scale_colour_manual(values = site_colours) +
#     scale_fill_manual(values = site_colours) +
#     geom_point(size = 4) +
#     geom_ribbon(data = site_preds,
#                 mapping = aes(x = year,
#                               ymin = pred_l95,
#                               ymax = pred_u95,
#                               fill = spp_phen,
#                               group = spp_phen),
#                 alpha = 0.5,
#                 inherit.aes = FALSE) +
#     geom_line(data = site_preds,
#               mapping = aes(x = year,
#                             y = pred,
#                             colour = spp_phen,
#                             group = spp_phen),
#               inherit.aes = FALSE,
#               size = 1) +
#     scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
#     scale_y_continuous(limits = c(90, 250), breaks = seq(100, 240, 20)) +
#     xlab(label = "") +
#     ylab(label = ylab_filter("Mean Phenology (DoY)", site_to_plot)) +
#     annotate("text", x = 2014, y = 245, label = site_to_plot, 
#              colour = site_colour, size = 7, hjust = 1) +
#     theme_bw() +
#     theme(panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           axis.text = element_text(colour = "black", size = 15),
#           axis.title = element_text(size = 15),
#           legend.position = "none") 
# }
# # plot all sites
# list2env(lapply(setNames(site_names, 
#                          make.names(paste0(site_names,
#                                            "_phen_plot"))),
#                 plot_phen), 
#          envir = .GlobalEnv)
# 
# # Set panel layout
# panel_layout <- rbind(c(1,2,3,4))
# 
# # Combine plots and export
# plot_grob <- grid.arrange(ALEXFIORD_phen_plot, 
#                           BARROW_phen_plot, 
#                           QHI_phen_plot, 
#                           ZACKENBERG_phen_plot,
#                           layout_matrix = panel_layout)
# ggsave(filename = paste0(script_path, "/coastal_phen_plot_MCMCglmm_5d.png"), 
#        plot = plot_grob, scale = 1.7, 
#        width = 10, 
#        height = (10/3), 
#        dpi = 600)
# 
# ### Calculate decadal advance ---- 
# phen_trends <- phen_trends %>% mutate(adv_per_dec = slope * 10)
# 
# # Min and max advance
# phen_trends %>% summarise(min_adv_per_dec = round(min(adv_per_dec),2), 
#                           max_adv_per_dec = round(max(adv_per_dec),2))
# 
# # Export
# write.csv(phen_trends, file = paste0(script_path, "coastal_phen_trends_5d.csv"), row.names = F)
# 
# 
# 
# 
# ####### Aditional calcs: 2) Begining of season data -20 instead of -10 ----
# 
# rm(list = ls())
# # Housekeeping
# script_path <- "scripts/users/jassmann/phenology/coastal_phenology/"
# load(file =paste0(script_path, "coastal_phen.Rda"))
# 
# # remove NAs (211 values)
# coastal_phen <- coastal_phen[!is.na(coastal_phen$doy),]
# 
# # Site names
# site_names <- unique(as.character(coastal_phen$site_name))
# 
# # Site spp phen combos
# site_spp_phen_combos <- unique(coastal_phen$site_spp_phen)
# 
# # Spp phen combos
# spp_phen_combos <- sort(unique(coastal_phen$site_spp_phen))
# 
# # Set colours
# colour_theme_sites <- data.frame(site_name = c(site_names),
#                                  colour = c("#324D5CFF", 
#                                             "#46B29DFF", 
#                                             "#C2A33EFF",
#                                             "#E37B40FF")) # 5th colour #F53855FF
# colour_theme_spp_phen <- data.frame(site_name = c(rep(site_names[1], 8), 
#                                                   rep(site_names[2], 7),
#                                                   rep(site_names[3], 3),
#                                                   rep(site_names[4], 6)),
#                                     site_spp_phen = spp_phen_combos,
#                                     spp_phen = c(unique(as.character(coastal_phen[coastal_phen$site_name == site_names[1],]$spp_phen)),
#                                                  unique(as.character(coastal_phen[coastal_phen$site_name == site_names[2],]$spp_phen)),
#                                                  unique(as.character(coastal_phen[coastal_phen$site_name == site_names[3],]$spp_phen)),
#                                                  unique(as.character(coastal_phen[coastal_phen$site_name == site_names[4],]$spp_phen))),
#                                     colour = c(
#                                       colorRampPalette(c("#19262EFF","#84CBF2FF"), 
#                                                        interpolate = "linear")(8), # ALEXFIIORD
#                                       colorRampPalette(c("#112B26FF","#5BE8CDFF"), 
#                                                        interpolate = "linear")(7), # BARROW
#                                       colorRampPalette(c("#382F12FF","#C2A33EFF"), 
#                                                        interpolate = "linear")(3), # QHI
#                                       colorRampPalette(c("#422413FF","#FF8A48FF"), 
#                                                        interpolate = "linear")(6) # ZACKENBERG
#                                     ), stringsAsFactors = F)
# 
# # Helper function that returns ylab only for Alexfiord (first site)
# ylab_filter <- function(label, site_to_plot){
#   if(site_to_plot == "ALEXFIORD") {
#     return(label)
#   } else {
#     return("")
#   }
# }
# 
# # Define function to run MCMCglmm models
# trend_model <- function(site_spp_sphen_select){
#   cphen_subset <- coastal_phen %>% 
#     filter(site_spp_phen == site_spp_sphen_select) 
#   model <- MCMCglmm(cbind(prior_visit, doy) ~ year, data = as.data.frame(cphen_subset),
#                     family = "cengaussian",
#                     nitt = 20000)
#   return(model)
# }
# 
# 
# # adjust data set
# # NB re-run begining of script to re-load data
# coastal_phen[(coastal_phen$doy - coastal_phen$prior_visit) == 10,
#              ]$prior_visit <- 
#   coastal_phen[(coastal_phen$doy - coastal_phen$prior_visit) == 10,
#                ]$doy - 20
# 
# # Execute over all site_spp_phen_combinations
# model_list <- lapply(site_spp_phen_combos, trend_model)
# names(model_list) <- make.names(paste0(site_spp_phen_combos, "_model"))
# list2env(model_list, envir = .GlobalEnv)                                
# 
# ### Gather outputs 
# ## Prep data frame
# phen_trends <- data.frame(site_spp_phen = site_spp_phen_combos,
#                           intercept = NA,
#                           intercept_l95 = NA,
#                           intercept_u95 = NA,          
#                           intercept_esmpl = NA,
#                           intercept_pMCMC = NA,
#                           slope = NA,
#                           slope_l95 = NA,
#                           slope_u95 = NA,
#                           slope_esmpl = NA,
#                           slope_pMCMC = NA)
# ## Define funciton to extract slope parameters and estimates 
# extr_outputs <- function(site_spp_sphen_select){
#   model <- get(paste0(site_spp_sphen_select, "_model"))
#   model_sum <- summary(model)
#   # retrive row index for output dataframe
#   row_no <- which(phen_trends$site_spp_phen == site_spp_sphen_select)
#   # use one liners to extract intersept and slope estimators
#   phen_trends[row_no,2:6] <<- model_sum$solutions[1,]
#   phen_trends[row_no,7:11] <<- model_sum$solutions[2,]
# }
# 
# ## Apply to all site_spp_phen combos
# lapply(site_spp_phen_combos, extr_outputs)
# 
# ### Claculate predictions
# # prep data frame
# phen_predicts <- data.frame(site_spp_phen = coastal_phen %>% 
#                               select(site_spp_phen, year) %>%
#                               distinct(site_spp_phen, year) %>%
#                               select(site_spp_phen) %>%
#                               unlist() %>% 
#                               as.character(),
#                             year = coastal_phen %>% 
#                               select(site_spp_phen, year) %>%
#                               distinct(site_spp_phen, year) %>%
#                               select(year) %>%
#                               unlist() %>% 
#                               as.numeric(),
#                             pred = NA,
#                             pred_l95 = NA,
#                             pred_u95 = NA)
# 
# # Define function to extract predicitons
# extract_pred <- function(site_spp_phen_select){
#   model <- get(paste0(site_spp_phen_select, "_model"))
#   row_nos <- which(phen_predicts$site_spp_phen == site_spp_phen_select)
#   preds <- predict.MCMCglmm(model, interval = "confidence", type = "response")
#   preds <- unique(data.frame(preds, year = as.numeric(model$X[,2])))
#   preds <- preds[,c(4,1,2,3)]
#   phen_predicts[row_nos, 3:5] <<- preds[,2:4]
#   paste0("Done with: ", site_spp_phen_select, " ;")
# }
# lapply(site_spp_phen_combos, extract_pred)
# 
# # create site name and spp_phen columns
# phen_predicts$site_name <- gsub("_.*", "", phen_predicts$site_spp_phen) 
# phen_predicts$spp_phen <- gsub("^[A-Z]*_", "", phen_predicts$site_spp_phen) 
# 
# ### Plot phenology data with trends ----
# plot_phen <- function(site_to_plot){
#   
#   site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
#   site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
#   site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
#     group_by(spp_phen, year) %>% summarise(phen_mean = round(mean(doy, na.rm = T)))
#   
#   site_preds <- phen_predicts %>%
#     filter(site_name == site_to_plot)
#   
#   ggplot(data = site_phen, aes(x = year, 
#                                y = phen_mean, 
#                                colour = spp_phen, 
#                                fill= spp_phen, 
#                                group = spp_phen)) +
#     scale_colour_manual(values = site_colours) +
#     scale_fill_manual(values = site_colours) +
#     geom_point(size = 4) +
#     geom_line(data = site_preds,
#               mapping = aes(x = year,
#                             y = pred,
#                             colour = spp_phen,
#                             group = spp_phen),
#               inherit.aes = FALSE,
#               size = 1) +
#     geom_ribbon(data = site_preds,
#                 mapping = aes(x = year,
#                               ymin = pred_l95,
#                               ymax = pred_u95,
#                               fill = spp_phen,
#                               group = spp_phen),
#                 alpha = 0.5,
#                 inherit.aes = FALSE) +
#     scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
#     scale_y_continuous(limits = c(90, 250), breaks = seq(100, 240, 20)) +
#     xlab(label = "") +
#     ylab(label = ylab_filter("Mean Phenology (DoY)", site_to_plot)) +
#     annotate("text", x = 2014, y = 245, label = site_to_plot, 
#              colour = site_colour, size = 7, hjust = 1) +
#     theme_bw() +
#     theme(panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.line = element_line(colour = "black"),
#           axis.text = element_text(colour = "black", size = 15),
#           axis.title = element_text(size = 15),
#           legend.position = "none") 
# }
# # plot all sites
# list2env(lapply(setNames(site_names, 
#                          make.names(paste0(site_names,
#                                            "_phen_plot"))),
#                 plot_phen), 
#          envir = .GlobalEnv)
# 
# # Set panel layout
# panel_layout <- rbind(c(1,2,3,4))
# 
# # Combine plots and export
# plot_grob <- grid.arrange(ALEXFIORD_phen_plot, 
#                           BARROW_phen_plot, 
#                           QHI_phen_plot, 
#                           ZACKENBERG_phen_plot,
#                           layout_matrix = panel_layout)
# ggsave(filename = paste0(script_path, "/coastal_phen_plot_MCMCglmm_20d.png"), 
#        plot = plot_grob, scale = 1.7, 
#        width = 10, 
#        height = (10/3), 
#        dpi = 600)
# 
# ### Calculate decadal advance ---- 
# phen_trends <- phen_trends %>% mutate(adv_per_dec = slope * 10)
# 
# # Min and max advance
# phen_trends %>% summarise(min_adv_per_dec = round(min(adv_per_dec),2), 
#                           max_adv_per_dec = round(max(adv_per_dec),2))
# 
# # Export
# write.csv(phen_trends, file = paste0(script_path, "coastal_phen_trends_20d.csv"), row.names = F)
# 
