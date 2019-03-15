# Visualise Choastal Phenology Data Trends
# Jakob Assmann j.assmann@ed.ac.uk 19 Feb 2018

# Dependencies
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(viridisLite)

## Housekeeping
script_path <- "scripts/users/jassmann/phenology/coastal_phenology/"
load(file = paste0(script_path, "coastal_phen.Rda"))
# load(file = paste0(script_path, "models/coastal_phen_rslopes_full_2018-03-21.Rda"))
load(file = "/Volumes/csce/biology/users/s1043792/scratch/coastal_phen_rslopes_full_run3_2018-03-21.Rda")
coastal_phen_rslopes_full <- best_run
rm(best_run)

## Set colour themes
# There are four sites
site_names <- unique(as.character(coastal_phen$site_name))
colour_theme_sites <- data.frame(site_name = c(site_names),
                                 name_pretty =  c("Alexandra Fiord", 
                                                  "Utqiaġvik",
                                                  "Qikiqtaruk",
                                                  "Zackenberg"),
                                 colour = c("#324D5CFF", 
                                            "#46B29DFF", 
                                            "#C2A33EFF",
                                            "#E37B40FF")) # 5th colour #F53855FF
spp_phen_combos <- sort(unique(coastal_phen$site_spp_phen))
# Set colours 
colour_theme_spp_phen <- data.frame(site_name = c(rep(site_names[1], 8), rep(site_names[2], 7), rep(site_names[3], 3), rep(site_names[4], 6)),
                                    site_spp_phen = spp_phen_combos,
                                    spp_phen = c(unique(as.character(coastal_phen[coastal_phen$site_name == site_names[1],]$spp_phen)),
                                                 unique(as.character(coastal_phen[coastal_phen$site_name == site_names[2],]$spp_phen)),
                                                 unique(as.character(coastal_phen[coastal_phen$site_name == site_names[3],]$spp_phen)),
                                                 unique(as.character(coastal_phen[coastal_phen$site_name == site_names[4],]$spp_phen))),
                                    colour = c(
                                      colorRampPalette(c("#19262EFF","#84CBF2FF"), interpolate = "linear")(8), # ALEXFIIORD
                                      colorRampPalette(c("#112B26FF","#5BE8CDFF"), interpolate = "linear")(7), # BARROW
                                      colorRampPalette(c("#382F12FF","#C2A33EFF"), interpolate = "linear")(3), # QHI
                                      colorRampPalette(c("#422413FF","#FF8A48FF"), interpolate = "linear")(6) # ZACKENBERG
                                    ), stringsAsFactors = F)

## Helper function that returns ylab only for Alexfiord (first site)
ylab_filter <- function(label, site_to_plot){
  if(site_to_plot == "ALEXFIORD") {
    return(label)
  } else {
    return("")
  }
}

### Preparation of posterior distributions for slope parameters of each site_spp_phen combination ----

# Prep data frame
post_dists <- data.frame(
  site_name = c(rep(site_names[1], 8), rep(site_names[2], 7), rep(site_names[3], 3), rep(site_names[4], 6)),
  site_spp_phen = spp_phen_combos,
  snowmelt_mean = mean(coastal_phen_rslopes_full$Sol[,"snowmelt_rel_scaled"]),
  snowmelt_mean_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"snowmelt_rel_scaled"], 0.05)),
  snowmelt_mean_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"snowmelt_rel_scaled"], 0.95)),
  phase_temp_mean = mean(coastal_phen_rslopes_full$Sol[,"phase_temp_rel_scaled"]), 
  phase_temp_mean_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"phase_temp_rel_scaled"], 0.05)),
  phase_temp_mean_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"phase_temp_rel_scaled"], 0.95)),
  onset_ice_melt_mean = mean(coastal_phen_rslopes_full$Sol[,"onset_ice_melt_rel_scaled"]),
  onset_ice_melt_mean_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"onset_ice_melt_rel_scaled"], 0.05)),
  onset_ice_melt_mean_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"onset_ice_melt_rel_scaled"], 0.95)),
  year_mean = mean(coastal_phen_rslopes_full$Sol[,"year_rel_scaled"]),
  year_mean_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"year_rel_scaled"], 0.05)),
  year_mean_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,"year_rel_scaled"], 0.95)),
  snowmelt_spp = NA,
  snowmelt_spp_l95 = NA,
  snowmelt_spp_u95 = NA,
  phase_temp_spp = NA,
  phase_temp_spp_l95 = NA,
  phase_temp_spp_u95 = NA,
  onset_ice_melt_spp = NA,
  onset_ice_melt_spp_l95 = NA,
  onset_ice_melt_spp_u95 = NA,
  year_spp = NA,
  year_spp_l95 = NA,
  year_spp_u95 = NA,
  snowmelt_total = NA,
  snowmelt_total_l95 = NA,
  snowmelt_total_u95 = NA,
  phase_temp_total = NA,
  phase_temp_total_l95 = NA,
  phase_temp_total_u95 = NA,
  onset_ice_melt_total = NA,
  onset_ice_melt_total_l95 = NA,
  onset_ice_melt_total_u95 = NA,
  year_total = NA,
  year_total_l95 = NA,
  year_total_u95 = NA)

### Fill site_spp_phen specific columns
extract_post_dists <- function(site_spp_phen){
  
  # Extract site_spp_phen specific distribution parameters
  snowmelt_spp = mean(coastal_phen_rslopes_full$Sol[,paste0("snowmelt_rel_scaled.site_spp_phen.", site_spp_phen)])
  snowmelt_spp_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("snowmelt_rel_scaled.site_spp_phen.", site_spp_phen)], 0.05))
  snowmelt_spp_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("snowmelt_rel_scaled.site_spp_phen.", site_spp_phen)], 0.95))
  
  phase_temp_spp = mean(coastal_phen_rslopes_full$Sol[,paste0("phase_temp_rel_scaled.site_spp_phen.", site_spp_phen)])
  phase_temp_spp_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("phase_temp_rel_scaled.site_spp_phen.", site_spp_phen)], 0.05))
  phase_temp_spp_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("phase_temp_rel_scaled.site_spp_phen.", site_spp_phen)], 0.95))
  
  onset_ice_melt_spp = mean(coastal_phen_rslopes_full$Sol[,paste0("onset_ice_melt_rel_scaled.site_spp_phen.", site_spp_phen)])
  onset_ice_melt_spp_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("onset_ice_melt_rel_scaled.site_spp_phen.", site_spp_phen)], 0.05))
  onset_ice_melt_spp_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("onset_ice_melt_rel_scaled.site_spp_phen.", site_spp_phen)], 0.95))
  
  year_spp = mean(coastal_phen_rslopes_full$Sol[,paste0("year_rel_scaled.site_spp_phen.", site_spp_phen)])
  year_spp_l95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("year_rel_scaled.site_spp_phen.", site_spp_phen)], 0.05))
  year_spp_u95 = as.numeric(quantile(coastal_phen_rslopes_full$Sol[,paste0("year_rel_scaled.site_spp_phen.", site_spp_phen)], 0.95))
  
  # Export to post_dists df
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$snowmelt_spp <<- snowmelt_spp
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$snowmelt_spp_l95 <<- snowmelt_spp_l95
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$snowmelt_spp_u95 <<- snowmelt_spp_u95
  
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$phase_temp_spp <<- phase_temp_spp
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$phase_temp_spp_l95 <<- phase_temp_spp_l95
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$phase_temp_spp_u95 <<- phase_temp_spp_u95
  
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$onset_ice_melt_spp <<- onset_ice_melt_spp
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$onset_ice_melt_spp_l95 <<- onset_ice_melt_spp_l95
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$onset_ice_melt_spp_u95 <<- onset_ice_melt_spp_u95
  
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$year_spp <<- year_spp
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$year_spp_l95 <<- year_spp_l95
  post_dists[post_dists$site_spp_phen == site_spp_phen,]$year_spp_u95 <<- year_spp_u95
  
  return(paste0("Finished ", site_spp_phen, "."))
}

# Execute function across all spp
sapply(post_dists$site_spp_phen, extract_post_dists)

### Calculate total effect sizes and credible intervals
post_dists$snowmelt_total = post_dists$snowmelt_mean + 
  post_dists$snowmelt_spp
post_dists$snowmelt_total_l95 = post_dists$snowmelt_mean_l95 + 
  post_dists$snowmelt_spp_l95
post_dists$snowmelt_total_u95 = post_dists$snowmelt_mean_u95 +
  post_dists$snowmelt_spp_u95

post_dists$phase_temp_total = post_dists$phase_temp_mean + 
  post_dists$phase_temp_spp
post_dists$phase_temp_total_l95 = post_dists$phase_temp_mean_l95 +
  post_dists$phase_temp_spp_l95
post_dists$phase_temp_total_u95 = post_dists$phase_temp_mean_u95 +
  post_dists$phase_temp_spp_u95

post_dists$onset_ice_melt_total = post_dists$onset_ice_melt_mean + 
  post_dists$onset_ice_melt_spp
post_dists$onset_ice_melt_total_l95 = post_dists$onset_ice_melt_mean_l95 +
  post_dists$onset_ice_melt_spp_l95
post_dists$onset_ice_melt_total_u95 = post_dists$onset_ice_melt_mean_u95 +
  post_dists$onset_ice_melt_spp_u95

post_dists$year_total = post_dists$year_mean +
  post_dists$year_spp
post_dists$year_total_l95 = post_dists$year_mean_l95 +
  post_dists$year_spp_l95
post_dists$year_total_u95 = post_dists$year_mean_u95 + 
  post_dists$year_spp_u95

# Turn into long form for plotting

# Subset and expand to long form
post_dists_total <- gather(post_dists[,c("site_name", 
                                         "site_spp_phen",
                                         "snowmelt_total",
                                         "phase_temp_total", 
                                         "onset_ice_melt_total")],
                           key = "predictor", value = "effect_size", 3:5)
post_dists_lower <- gather(post_dists[,c("site_name", 
                                         "site_spp_phen",
                                         "snowmelt_total_l95", 
                                         "phase_temp_total_l95", 
                                         "onset_ice_melt_total_l95")],
                           key = "predictor", value = "lower", 3:5)
post_dists_upper <- gather(post_dists[,c("site_name", 
                                         "site_spp_phen",
                                         "snowmelt_total_u95",
                                         "phase_temp_total_u95", 
                                         "onset_ice_melt_total_u95")], 
                           key = "predictor", value = "upper", 3:5)
# Combine dfs and clear up and rename
post_dists_total <- cbind(post_dists_total, 
                          post_dists_lower[,"lower"],
                          post_dists_upper[,"upper"])
rm(post_dists_lower, post_dists_upper)
names(post_dists_total) <- c("site_name", 
                             "site_spp_phen",
                             "predictor",
                             "effect_size",
                             "lower", 
                             "upper")
post_dists_total[post_dists_total$predictor == "snowmelt_total",
                 ]$predictor <- "Snowmelt"
post_dists_total[post_dists_total$predictor == "phase_temp_total",
                 ]$predictor <- "Temperature"
post_dists_total[post_dists_total$predictor == "onset_ice_melt_total",
                 ]$predictor <- "Sea-Ice"

# Create extra dummy species so that plotting is equal. 
# NB ggplot sorts as levels are specified in factors,
# so I use underscores in the dummy species to manipulate the order the non dummy species are displayed in
dummy_spp <- data.frame(site_name = rep(c(rep("ALEXFIORD", 2),
                                          rep("BARROW",3),
                                          rep("QHI", 7), 
                                          rep("ZACKENBERG", 4)), 3),
                        site_spp_phen = rep(c("_ALEXFIORD1",
                                              "ALEXFIORD2",
                                              "_BARROW1",
                                              "_BARROW2",
                                              "BARROW3", 
                                              "_QHI1", 
                                              "_QHI2",
                                              "_QHI3", 
                                              "_QHI4", "QHI5",
                                              "QHI6",
                                              "QHI7", 
                                              "_ZACKENBERG1", 
                                              "_ZACKENBERG2",
                                              "ZACKENBERG3", 
                                              "ZACKENBERG4"),3),
                        predictor = c(rep("Snowmelt", 16), 
                                      rep("Temperature", 16), 
                                      rep("Sea-Ice",16)),
                        effect_size = 0,
                        lower = 0,
                        upper = 0)

post_dists_total <- rbind(post_dists_total, dummy_spp)
# Refesh order in factor, first "unfactor"
post_dists_total$site_spp_phen <- as.character(post_dists_total$site_spp_phen)
# Sort df aplphabetically
post_dists_total <- post_dists_total[order(post_dists_total$site_spp_phen),]
# Re-factor
post_dists_total$site_spp_phen <- factor(post_dists_total$site_spp_phen)
# Re-order predictor factor for plottin (snomwlet, temp then sea ice)
post_dists_total$predictor <- factor(post_dists_total$predictor)
post_dists_total$predictor <- factor(post_dists_total$predictor, 
                                     levels = levels(post_dists_total$predictor)[c(2,3,1)])
post_dists_total$predictor

# Add transparent colours to the colour frame
dummy_spp <- data.frame(site_name = c(rep("ALEXFIORD", 2), 
                                      rep("BARROW",3), 
                                      rep("QHI", 7), 
                                      rep("ZACKENBERG", 4)),
                        site_spp_phen = c("_ALEXFIORD1",
                                          "ALEXFIORD2",
                                          "_BARROW1", 
                                          "_BARROW2", 
                                          "BARROW3",
                                          "_QHI1",
                                          "_QHI2", 
                                          "_QHI3", 
                                          "_QHI4",
                                          "QHI5", 
                                          "QHI6",
                                          "QHI7", 
                                          "_ZACKENBERG1",
                                          "_ZACKENBERG2",
                                          "ZACKENBERG3",
                                          "ZACKENBERG4"),
                        spp_phen = c(paste0("ALEXFIORD", 
                                            as.character(seq(1,2))),
                                     paste0("BARROW", as.character(seq(1,3))), 
                                     paste0("QHI", as.character(seq(1,7))), 
                                     paste0("ZACKENBERG", as.character(seq(1,4)))),
                        colour = "#FFFFFF00")
colour_theme_spp_phen <- rbind(colour_theme_spp_phen, dummy_spp)
# Bring in same order as post_dists df
colour_theme_spp_phen$site_spp_phen <- as.character(colour_theme_spp_phen$site_spp_phen)
colour_theme_spp_phen <- colour_theme_spp_phen[order(colour_theme_spp_phen$site_spp_phen),]
colour_theme_spp_phen$site_spp_phen <- factor(colour_theme_spp_phen$site_spp_phen)


plot_spp_effect_sizes <- function(site_to_plot){
  
  site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  name_pretty <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$name_pretty
  spp_dists <- post_dists_total %>% filter(site_name == site_to_plot)
  
  ggplot(data = spp_dists, aes(x = as.integer(predictor),
                               y = effect_size, 
                               ymin = lower, 
                               ymax = upper, 
                               colour = site_spp_phen,
                               group = site_spp_phen)) +
    scale_colour_manual(values = site_colours) +
    geom_point(stat="identity", position = position_dodge(width = 1), size = 3) +
    geom_errorbar(position = position_dodge(width = 1), width = 1) +
    scale_x_continuous(labels=c("Snowmelt","Temperature", "Sea-Ice"),
                       breaks=c(1,2,3), limits=c(0.5,3.5), expand=c(0,0)) +
    scale_y_continuous(limits = c(-8, 8), breaks = seq(-8,8, 2)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Effect Size (Scaled)", site_to_plot)) +
    annotate("text", x = 3.4, y = 8, label = name_pretty, 
             colour = site_colour, size = 7, hjust = 1) +
    annotate("segment", x = 0.5, xend = 3.5, y = 0, yend = 0) +
    annotate("segment", x = 1.5, xend = 1.5, y = -8, yend = 7, linetype = 2) +
    annotate("segment", x = 2.5, xend = 2.5, y = -8, yend = 7, linetype = 2) +
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
                                           "_spp_effect_size_plot"))),
                plot_spp_effect_sizes), 
         envir = .GlobalEnv)
# plot
panel_layout <- rbind(c(1,2,3,4))

# png(filename = paste0(script_path, "/coastal_spp_effect_size_plot.png"), width = 4800, height = 1440, res = 288)
# grid.arrange(ALEXFIORD_spp_effect_size_plot, 
#              BARROW_spp_effect_size_plot, 
#              QHI_spp_effect_size_plot, 
#              ZACKENBERG_spp_effect_size_plot,
#              layout_matrix = panel_layout)
# dev.off()

# Now as EPS
plot_grob <- grid.arrange(ALEXFIORD_spp_effect_size_plot, 
                          BARROW_spp_effect_size_plot, 
                          QHI_spp_effect_size_plot, 
                          ZACKENBERG_spp_effect_size_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_spp_effect_size_plot.png"), 
       plot = plot_grob, scale = 1.7, 
       width = 10, 
       height = (10/3), 
       dpi = 600)
# clean up
rm(list=ls()[grepl('*_plot*',ls())])

### Violin plots of mean effect sizes and distributons
# Extract data summary (will speed up later data extraction)
effect_sizes <- data.frame(predictor = c("Snowmelt", "Temperature", "Sea-Ice"),
                           pred_mean = c(summary(coastal_phen_rslopes_full)$solutions["snowmelt_rel_scaled","post.mean"],
                                         summary(coastal_phen_rslopes_full)$solutions["phase_temp_rel_scaled","post.mean"],
                                         summary(coastal_phen_rslopes_full)$solutions["onset_ice_melt_rel_scaled","post.mean"]),
                           pred_var = c(summary(coastal_phen_rslopes_full)$Gcovariances["snowmelt_rel_scaled:snowmelt_rel_scaled.site_spp_phen","post.mean"],
                                         summary(coastal_phen_rslopes_full)$Gcovariances["phase_temp_rel_scaled:phase_temp_rel_scaled.site_spp_phen","post.mean"],
                                         summary(coastal_phen_rslopes_full)$Gcovariances["onset_ice_melt_rel_scaled:onset_ice_melt_rel_scaled.site_spp_phen","post.mean"])
                           )
# Sample random normal distirbution with parameters
# To fill violin plot
violin_data <- data.frame(snowmelt_sample = rnorm(n = 10000,
                                                  mean = effect_sizes[effect_sizes$predictor == "Snowmelt",]$pred_mean,
                                                  sd = sqrt(effect_sizes[effect_sizes$predictor == "Snowmelt",]$pred_var)),
                          temp_sample = rnorm(n = 10000,
                                                  mean = effect_sizes[effect_sizes$predictor == "Temperature",]$pred_mean,
                                                  sd = sqrt(effect_sizes[effect_sizes$predictor == "Temperature",]$pred_var)),
                          ice_sample = rnorm(n = 10000,
                                                  mean = effect_sizes[effect_sizes$predictor == "Sea-Ice",]$pred_mean,
                                                  sd = sqrt(effect_sizes[effect_sizes$predictor == "Sea-Ice",]$pred_var)))
# gather into long form
violin_data <- gather(violin_data, key = predictor, value = effect_size, 1:3)
violin_data[violin_data$predictor == "snowmelt_sample",]$predictor <- "Snowmelt"
violin_data[violin_data$predictor == "temp_sample",]$predictor <- "Temperature"
violin_data[violin_data$predictor == "ice_sample",]$predictor <- "Sea-Ice"
violin_data$predictor <- factor(violin_data$predictor)
violin_data$predictor <- factor(violin_data$predictor, levels = levels(factor(violin_data$predictor))[c(2,3,1)])

violin_plot <- ggplot(data = violin_data, aes(x = predictor, y = effect_size, fill = predictor)) +
  geom_violin() +
  scale_fill_manual(values = viridis(3)[c(2,3,1)]) +
  #scale_x_continuous(labels=c("Snowmelt","Temperature", "Sea-Ice"), breaks=c(1,2,3), limits=c(0.2,3.8), expand=c(0,0)) +
  scale_y_continuous(limits = c(-7, 10), breaks = seq(-8,10, 2)) +
  annotate("segment", x = 0.75, xend = 1.25, 
           y = effect_sizes[effect_sizes$predictor == "Snowmelt",]$pred_mean, 
           yend = effect_sizes[effect_sizes$predictor == "Snowmelt",]$pred_mean, 
           size = 1.2) +
  annotate("segment", x = 1.72, xend = 2.28,
           y = effect_sizes[effect_sizes$predictor == "Temperature",]$pred_mean,
           yend = effect_sizes[effect_sizes$predictor == "Temperature",]$pred_mean,
           size = 1.2) +
  annotate("segment", x = 2.55, xend = 3.45,
           y = effect_sizes[effect_sizes$predictor == "Sea-Ice",]$pred_mean, 
           yend = effect_sizes[effect_sizes$predictor == "Sea-Ice",]$pred_mean, 
           size = 1.2) +
  annotate("segment", x = 0, xend = 4, y = 0, yend = 0, linetype = 2) +
  xlab(label = "") +
  ylab(label = "Effect Size (Scaled)") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(size = 20),
        legend.position = "none")

ggsave(paste0(script_path, "coastal_phen_violin.png"), violin_plot)

# Export effect sizes in tabular format
post_dists_export <- post_dists_total %>% filter(effect_size != 0)
write.csv(post_dists_export[order(post_dists_export$predictor),], 
          file = paste0(script_path, "post_dists.csv"), row.names = F)

# Export table of all site spp phen combinations with associated colours
colour_export <- colour_theme_spp_phen[colour_theme_spp_phen$colour != "#FFFFFF00",-2]
itex.data <- read.csv("scripts/users/jassmann/phenology/phenology_data/ITEX_data_Janet/CCIN12722_20171116_Arctic_phenology_database_1992-2014.csv")
colour_export$species <- sapply(substr(colour_export$spp_phen, 1, 6),
                                function(x) {as.character(first(itex.data[itex.data$spp == x,]$species))})
colour_export$genus <- sapply(substr(colour_export$spp_phen, 1, 6),
                              function(x) {as.character(first(itex.data[itex.data$spp == x,]$genus))})
colour_export$phen_stage <- substring(colour_export$spp_phen, 8)
colour_export <- colour_export[,c(1,4,5,6,3)]
write.csv(colour_export,
          file = paste0(script_path, "site_spp_phen_colour.csv"),
          row.names = F)

# Make a pretty legend plot
colour_export[colour_export$phen_stage == "green_up",]$phen_stage <- "green up"
legend_plot <- function(site_to_plot){
  site_colours <- colour_export[colour_export$site_name == site_to_plot,]
  site_colours$spp <- paste0(substr(site_colours$genus, 1,1), ". ", site_colours$species,
                             " ", site_colours$phen_stage)
  dummy_data <- data.frame(year = seq(1990,2016),
                           values = rnorm(27) * 70 + 170) # Create some dummy data actually not needed.
  
  ggplot(data = dummy_data, aes(x = year,
                               y = values)) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5), expand = c(0,0)) +
    scale_y_continuous(limits = c(170, 250), breaks = seq(100, 240, 20), expand = c(0,0)) +
    xlab(label = "") +
    ylab(label = "") +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 244, ymax = 250, fill = site_colours$colour[1],
             colour = site_colours$colour[1]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 234, ymax = 240, fill = site_colours$colour[2],
             colour = site_colours$colour[2]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 224, ymax = 230, fill = site_colours$colour[3],
             colour = site_colours$colour[3]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 214, ymax = 220, fill = site_colours$colour[4],
             colour = site_colours$colour[4]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 204, ymax = 210, fill = site_colours$colour[5],
             colour = site_colours$colour[5]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 194, ymax = 200, fill = site_colours$colour[6],
             colour = site_colours$colour[6]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 184, ymax = 190, fill = site_colours$colour[7],
             colour = site_colours$colour[7]) +
    annotate("rect", xmin = 1990, xmax = 1992,
             ymin = 174, ymax = 180, fill = site_colours$colour[8],
             colour = site_colours$colour[8]) +
    
    annotate("text", x = 1993, y = 244, label = site_colours$spp[1], 
             size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 234, label = site_colours$spp[2], 
             size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 224, label = site_colours$spp[3], 
             size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 214, label = site_colours$spp[4],
             size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 204, label = site_colours$spp[5],
              size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 194, label = site_colours$spp[6], 
            size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 184, label = site_colours$spp[7],
             size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    annotate("text", x = 1993, y = 174, label = site_colours$spp[8],
             size = 6, hjust = 0, vjust = 0, fontface = "italic") +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "white"),
          axis.text = element_text(colour = "white", size = 15),
          axis.title = element_text(size = 20),
          axis.ticks = element_blank(),
          legend.position = "none")
}

list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_legend_plot"))),
                legend_plot), 
         envir = .GlobalEnv)
plot_grob <- grid.arrange(ALEXFIORD_legend_plot, 
                          BARROW_legend_plot, 
                          QHI_legend_plot, 
                          ZACKENBERG_legend_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_legend.png"),
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3)/2, 
       dpi = 600)
