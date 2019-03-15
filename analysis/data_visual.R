# Visualise Choastal Phenology Data Trends
# Jakob Assmann j.assmann@ed.ac.uk 19 Feb 2018

# Dependencies
library(dplyr)
library(ggplot2)
library(gridExtra)

# Housekeeping
script_path <- "scripts/users/jassmann/phenology/coastal_phenology/"
load(file =paste0(script_path, "coastal_phen.Rda"))

# There are four sites
site_names <- unique(as.character(coastal_phen$site_name))
colour_theme_sites <- data.frame(site_name = c(site_names),
                          colour = c("#324D5CFF", 
                                     "#46B29DFF", 
                                     "#C2A33EFF",
                                     "#E37B40FF")) # 5th colour #F53855FF
# qhi original F0CA4D

# There are 24 species - phenology phase combinations in the dataset
# ALEXFIORD has 8
# BARROW has 7
# QHI has 3
# Zackenberg has 6
coastal_phen %>% group_by(site_name) %>% 
  summarise(n_comobos = length(unique(site_spp_phen)))
spp_phen_combos <- sort(unique(coastal_phen$site_spp_phen))
# Set colours 
colour_theme_spp_phen <- data.frame(site_name = c(rep(site_names[1], 8), 
                                                  rep(site_names[2], 7),
                                                  rep(site_names[3], 3),
                                                  rep(site_names[4], 6)),
                                    site_spp_phen = spp_phen_combos,
                                    spp_phen = c(unique(as.character(coastal_phen[coastal_phen$site_name == site_names[1],]$spp_phen)),
                                                 unique(as.character(coastal_phen[coastal_phen$site_name == site_names[2],]$spp_phen)),
                                                 unique(as.character(coastal_phen[coastal_phen$site_name == site_names[3],]$spp_phen)),
                                                 unique(as.character(coastal_phen[coastal_phen$site_name == site_names[4],]$spp_phen))),
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
## Visualise temperature 
# ----
# phen_phase temperature plot function
plot_phase_temp <- function(site_to_plot){
  site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(site_spp_phen, year) %>% summarise(phase_temp = first(phase_temp))
  
  ggplot(data = site_phen, aes(x = year, 
                               y = phase_temp, 
                               colour = site_spp_phen, 
                               group = site_spp_phen)) +
    scale_colour_manual(values = site_colours) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = F) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(-10, 5), breaks = seq(-10,5, 5)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Phase Temperature (C)", site_to_plot)) +
    annotate("text", x = 2014, y = 4.7, label = site_to_plot, 
             colour = site_colour, size = 7, hjust = 1) +
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
                                           "_phase_temp_plot"))),
                plot_phase_temp), 
         envir = .GlobalEnv)
# plot
panel_layout <- rbind(c(1,2,3,4))
# png(filename = paste0(script_path, "/coastal_phase_temp_plot.png"), width = 2400, height = 720)
# grid.arrange(ALEXFIORD_phase_temp_plot, 
#              BARROW_phase_temp_plot, 
#              QHI_phase_temp_plot, 
#              ZACKENBERG_phase_temp_plot,
#              layout_matrix = panel_layout)
# dev.off()
plot_grob <- grid.arrange(ALEXFIORD_phase_temp_plot, 
                          BARROW_phase_temp_plot, 
                          QHI_phase_temp_plot, 
                          ZACKENBERG_phase_temp_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_phase_temp_plot.png"), 
       plot = plot_grob,
       scale = 1.7,
       width = 10, 
       height = (10/3),
       dpi = 600)


# clean up
rm(list=ls()[grepl('*_plot*',ls())])

# site_phen_phase temperature plot function
plot_site_phase_temp <- function(site_to_plot){
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(year) %>% summarise(spring_temp = first(site_phase_temp))
  ggplot(data = site_phen, aes(x = year, y = spring_temp)) +
    geom_point(size = 4, colour = site_colour) +
    geom_smooth(method = "lm", colour = site_colour) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(-10, 5), breaks = seq(-10,5, 5)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Spring Temperature (C)", site_to_plot)) +
    annotate("text", x = 2014, y = 4.7, label = site_to_plot, 
             colour = site_colour, size = 7, hjust = 1) +
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
                                           "_site_phase_temp_plot"))),
                plot_site_phase_temp), 
         envir = .GlobalEnv)
# plot
panel_layout <- rbind(c(1,2,3,4))

# png(filename = paste0(script_path, "/coastal_site_phase_temp_plot.png"), width = 4800, height = 1440, res = 288)
# grid.arrange(ALEXFIORD_site_phase_temp_plot, 
#              BARROW_site_phase_temp_plot, 
#              QHI_site_phase_temp_plot, 
#              ZACKENBERG_site_phase_temp_plot,
#              layout_matrix = panel_layout)
# dev.off()

plot_grob <- grid.arrange(ALEXFIORD_site_phase_temp_plot, 
             BARROW_site_phase_temp_plot, 
             QHI_site_phase_temp_plot, 
             ZACKENBERG_site_phase_temp_plot,
             layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_site_phase_temp_plot.png"), 
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)

# clean up
rm(list=ls()[grepl('*_plot*',ls())])



# spring (April-June Average) temperature plot function
plot_spring_temp <- function(site_to_plot){
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(year) %>% summarise(spring_temp = first(spring_temp))
  ggplot(data = site_phen, aes(x = year, y = spring_temp)) +
    geom_point(size = 4, colour = site_colour) +
    geom_smooth(method = "lm", colour = site_colour) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(-10, 5), breaks = seq(-10,5, 5)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Spring Temperature (C)", site_to_plot)) +
    annotate("text", x = 2014, y = 4.7, label = site_to_plot, 
             colour = site_colour, size = 7, hjust = 1) +
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
                                           "_spring_temp_plot"))),
                plot_spring_temp), 
         envir = .GlobalEnv)
# plot
panel_layout <- rbind(c(1,2,3,4))

# png(filename = paste0(script_path, "/coastal_spring_temp_plot.png"), width = 4800, height = 1440, res = 288)
# grid.arrange(ALEXFIORD_spring_temp_plot, 
#              BARROW_spring_temp_plot, 
#              QHI_spring_temp_plot, 
#              ZACKENBERG_spring_temp_plot,
#              layout_matrix = panel_layout)
# dev.off()

plot_grob <- grid.arrange(ALEXFIORD_spring_temp_plot, 
                          BARROW_spring_temp_plot, 
                          QHI_spring_temp_plot, 
                          ZACKENBERG_spring_temp_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_spring_temp_plot.png"),
       plot = plot_grob,
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)

# clean up
rm(list=ls()[grepl('*_plot*',ls())])



## Visualise snowmelt 
# ----
# spring snowmelt plot function
plot_snowmelt <- function(site_to_plot){
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(year) %>% summarise(snowmelt = round(mean(snowmelt, na.rm = T)))
  ggplot(data = site_phen, aes(x = year, y = snowmelt)) +
    geom_point(size = 4, colour = site_colour) +
    geom_smooth(method = "lm", colour = site_colour) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(90, 250), breaks = seq(100, 240, 20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Snowmelt Date (DoY)", site_to_plot)) +
    annotate("text", x = 2014, y = 245, label = site_to_plot, 
             colour = site_colour, size = 7, hjust = 1) +
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

# png(filename = paste0(script_path, "/coastal_snowmelt_plot.png"), width = 4800, height = 1440, res = 288)
# grid.arrange(ALEXFIORD_snowmelt_plot, 
#              BARROW_snowmelt_plot, 
#              QHI_snowmelt_plot, 
#              ZACKENBERG_snowmelt_plot,
#              layout_matrix = panel_layout)
# dev.off()

plot_grob <- grid.arrange(ALEXFIORD_snowmelt_plot, 
                         BARROW_snowmelt_plot, 
                         QHI_snowmelt_plot, 
                         ZACKENBERG_snowmelt_plot,
                         layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_snowmelt_plot.png"),
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)


# clean up
rm(list=ls()[grepl('*_plot*',ls())])


# ----

## Visualise onset of sea-ice melt
# ----
# sea ice  plot function
site_to_plot <- site_names[1]
plot_sea_ice <- function(site_to_plot){
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(year) %>% summarise(onset_ice_melt = first(onset_ice_melt))
  ggplot(data = site_phen, aes(x = year, y = onset_ice_melt)) +
    geom_point(size = 4, colour = site_colour) +
    geom_smooth(method = "lm", colour = site_colour) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(90, 250), breaks = seq(100,240,20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Onset Sea-Ice Melt (DoY)", site_to_plot)) +
    annotate("text", x = 2014, y = 245, label = site_to_plot, colour = site_colour, size = 7, hjust = 1) +
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

# png(filename = paste0(script_path, "/coastal_sea_ice_plot.png"), width = 4800, height = 1440, res = 288)
# grid.arrange(ALEXFIORD_sea_ice_plot, 
#              BARROW_sea_ice_plot, 
#              QHI_sea_ice_plot, 
#              ZACKENBERG_sea_ice_plot,
#              layout_matrix = panel_layout)
# dev.off()

plot_grob <- grid.arrange(ALEXFIORD_sea_ice_plot, 
                          BARROW_sea_ice_plot, 
                          QHI_sea_ice_plot, 
                          ZACKENBERG_sea_ice_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_sea_ice_plot.png"),
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3),
       dpi = 600)


# clean up
rm(list=ls()[grepl('*_plot*',ls())])


## Phenoloy visualisation
# phen plot function
plot_phen <- function(site_to_plot){
  site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]$colour
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(spp_phen, year) %>% summarise(phen_mean = round(mean(doy, na.rm = T)))
  
  ggplot(data = site_phen, aes(x = year, 
                               y = phen_mean, 
                               colour = spp_phen, 
                               fill= spp_phen, 
                               group = spp_phen)) +
    scale_colour_manual(values = site_colours) +
    scale_fill_manual(values = site_colours) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = T, alpha = 0.25) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(90, 250), breaks = seq(100, 240, 20)) +
    xlab(label = "") +
    ylab(label = ylab_filter("Mean Phenology (DoY)", site_to_plot)) +
    annotate("text", x = 2014, y = 245, label = site_to_plot, 
             colour = site_colour, size = 7, hjust = 1) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black", size = 15),
          axis.title = element_text(size = 15),
          legend.position = "none")
}

legend_plot <- function(site_to_plot){
  site_colours <- colour_theme_spp_phen[colour_theme_spp_phen$site_name == site_to_plot,]
  site_colour <- colour_theme_sites[colour_theme_sites$site_name == site_to_plot,]$colour
  site_phen <- coastal_phen %>% filter(site_name == site_to_plot) %>% 
    group_by(spp_phen, year) %>% summarise(phen_mean = round(mean(doy, na.rm = T)))
  
  ggplot(data = site_phen, aes(x = year,
                               y = phen_mean,
                               colour = spp_phen,
                               group = spp_phen)) +
    scale_x_continuous(limits= c(1990, 2016), breaks = seq(1990, 2016, 5)) +
    scale_y_continuous(limits = c(90, 250), breaks = seq(100, 240, 20)) +
    xlab(label = "") +
    ylab(label = " ") +
    annotate("text", x = 1995, y = 240, label = site_colours$spp_phen[1], 
             colour = site_colours$colour[1], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 230, label = site_colours$spp_phen[2], 
             colour = site_colours$colour[2], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 220, label = site_colours$spp_phen[3], 
             colour = site_colours$colour[3], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 210, label = site_colours$spp_phen[4],
             colour = site_colours$colour[4], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 200, label = site_colours$spp_phen[5],
             colour = site_colours$colour[5], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 190, label = site_colours$spp_phen[6], 
             colour = site_colours$colour[6], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 180, label = site_colours$spp_phen[7],
             colour = site_colours$colour[7], size = 6, hjust = 0) +
    annotate("text", x = 1995, y = 170, label = site_colours$spp_phen[8],
             colour = site_colours$colour[8], size = 6, hjust = 0) +
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
# apply to all sites
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_phen_plot"))),
                plot_phen), 
         envir = .GlobalEnv)
list2env(lapply(setNames(site_names, 
                         make.names(paste0(site_names,
                                           "_legend_plot"))),
                legend_plot), 
         envir = .GlobalEnv)
# plot no legend
panel_layout <- rbind(c(1,2,3,4))

# png(filename = paste0(script_path, "/coastal_phen_plot.png"), width = 4800, height = 1440, res = 288)
# grid.arrange(ALEXFIORD_phen_plot, 
#              BARROW_phen_plot, 
#              QHI_phen_plot, 
#              ZACKENBERG_phen_plot,
#              layout_matrix = panel_layout)
# dev.off()

plot_grob <- grid.arrange(ALEXFIORD_phen_plot, 
                          BARROW_phen_plot, 
                          QHI_phen_plot, 
                          ZACKENBERG_phen_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_phen_plot.png"),
       plot = plot_grob, 
       scale = 1.7, 
       width = 10, 
       height = (10/3), 
       dpi = 600)


# plot with legend
panel_layout <- rbind(c(1,2,3,4),
                      c(5,6,7,8))

# png(filename = paste0(script_path, "/coastal_phen_plot_legend.png"), 
#     width = 1200, height = 720)
# grid.arrange(ALEXFIORD_phen_plot, 
#              BARROW_phen_plot, 
#              QHI_phen_plot, 
#              ZACKENBERG_phen_plot,
#              ALEXFIORD_legend_plot, 
#              BARROW_legend_plot, 
#              QHI_legend_plot, 
#              ZACKENBERG_legend_plot,
#              layout_matrix = panel_layout)
# dev.off()

plot_grob <- grid.arrange(ALEXFIORD_phen_plot, 
                          BARROW_phen_plot, 
                          QHI_phen_plot, 
                          ZACKENBERG_phen_plot,
                          ALEXFIORD_legend_plot, 
                          BARROW_legend_plot, 
                          QHI_legend_plot, 
                          ZACKENBERG_legend_plot,
                          layout_matrix = panel_layout)
ggsave(filename = paste0(script_path, "/coastal_phen_plot_legend.png"), 
       plot = plot_grob,
       scale = 1.7, 
       width = 10, 
       height = 2*(10/3),
       dpi = 600)


# clean up
rm(list=ls()[grepl('*_plot*',ls())])

