# Coastal phenology plant species and model fits.

# dependencies
library(ggplot2)
library(dplyr)
library(MCMCglmm)
library(gridExtra)
library(viridisLite)

# set folder path
script_path <- "scripts/users/jassmann/phenology/sea_ice_manuscript/janets_data/"

# load data
itex_data <- read.csv("scripts/users/jassmann/phenology/sea_ice_manuscript/janets_data/CCIN12722_20161221_Arctic_phenology_database_1992-2014.csv")

# data frame of species / site combos
spp_site_combos <- read.csv("scripts/users/jassmann/phenology/sea_ice_manuscript/janets_data/site_spp_combos.csv")

# subset data
phen_data <- itex_data %>% 
  inner_join(spp_site_combos[,c(3,4,7)], by = NULL) %>%
  group_by(site_name, spp, phenophase)

# mean phenology dates
summarise(phen_data, mean_phen = mean(day))

## wack out plots of mean annual phenology
phen_data_year_means <- phen_data %>% 
  ungroup() %>% 
  group_by(site_name, spp, phenophase, year) %>% 
  summarise(mean_phen = mean(day), mean_sd = sd(day)) %>%
  mutate(site_spp = paste0(site_name, "_", spp), 
         site_spp_pheno = paste0(site_name, "_", spp, "_", phenophase),
         spp_phen = paste0(spp, "_", phenophase))

# create plot funciton
site_plot <- function(selected_site){
  selected_site <- as.character(selected_site)
  temp_df <- phen_data_year_means %>% filter(site_name == selected_site)
  plot_colours <- viridis(5)[1:3]
  temp_plot <- ggplot(temp_df, aes(x = year, y = mean_phen, group = site_spp_pheno, colour = site_spp_pheno)) +
    geom_point() +
    geom_smooth(method = "lm") +
    scale_y_continuous(limits = c(140, 220), breaks = c(140, 160, 180, 200, 220)) +
    scale_x_continuous(limits = c(1995, 2016), breaks = c(1995, 2000, 2005, 2010, 2015)) +
    scale_colour_manual(values = plot_colours) +
    ylab("") +
    theme_bw() +
    annotate("text", x = 1995, y = 220, label = unique(temp_df$site_name), colour = "black", size = 5, hjust = 0) +
    annotate("text", x = 2010, y = 215, label = unique(temp_df$spp_phen)[1], colour = plot_colours[1], size = 4, hjust = 0.5 ) +
    annotate("text", x = 2010, y = 210, label = unique(temp_df$spp_phen)[2], colour = plot_colours[2], size = 4, hjust = 0.5 ) +
    annotate("text", x = 2010, y = 205, label = unique(temp_df$spp_phen)[3], colour = plot_colours[3], size = 4, hjust = 0.5 ) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(),
          axis.title = element_text(),
          legend.position = "none")
  return(temp_plot)
}

# plot all site spp combinations by site
list2env(lapply(setNames(unique(as.character(spp_site_combos$site_name)), 
                         make.names(paste0(unique(as.character(spp_site_combos$site_name)),
                                    "_phen_plot"))),
                site_plot), 
                envir = .GlobalEnv)

# make panel figure
panel_layout <- rbind(c(1,2,3,4))
png(filename = paste0(script_path, "/summary_plot.png"), width = 1000, height = 300)
grid.arrange(ALEXFIORD_phen_plot, BARROW_phen_plot, QHI_phen_plot, ZACKENBERG_phen_plot, layout_matrix = panel_layout)
dev.off()





####### Random Appendix Section
# Comparsion with QHI Data as there is a weird anomaly in Janets QHI data
library(tidyr)
qiki_phen_with_obs <- read.csv("scripts/phenology_scripts/data_2017/qiki_phen_with_before_2017.csv")
qhi_janet <- itex_data %>% filter(site_name == "QHI")
qiki_phen_with_obs_means <- qiki_phen_with_obs %>% 
  select(Spp, Year, 4:10) %>%
  group_by(Spp, Year) %>%
  summarise(p1_mean = mean(P1, na.rm = T),
            p2_mean = mean(P2, na.rm = T),
            p3_mean = mean(P3, na.rm = T),
            p4_mean = mean(P4, na.rm = T),
            p5_mean = mean(P5, na.rm = T),
            p6_mean = mean(P6, na.rm = T),
            p7_mean = mean(P7, na.rm = T))
qiki_means_long <- gather(qiki_phen_with_obs_means, "phen_stage", "doy", 3:9) %>% mutate(spp_phen = paste0(Spp, "_", phen_stage))            
SALARC_means_plot <- ggplot(qiki_means_long[qiki_means_long$Spp == "SALARC",], aes(x = Year, y = doy, group= spp_phen, colour =spp_phen)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_continuous(limits = c(140,250), breaks = c(140,160,180,200,220,240)) +
  scale_x_continuous(limits = c(2001,2017), breaks = c(2005,2010,2015)) +
  theme_bw()
ERIVAG_means_plot <- ggplot(qiki_means_long[qiki_means_long$Spp == "ERIVAG",], aes(x = Year, y = doy, group= spp_phen, colour =spp_phen)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_continuous(limits = c(130,170), breaks = c(130,140,150,160, 170)) +
  scale_x_continuous(limits = c(2001,2017), breaks = c(2005,2010,2015)) +
  theme_bw()
DRYINT_means_plot <- ggplot(qiki_means_long[qiki_means_long$Spp == "DRYINT",], aes(x = Year, y = doy, group= spp_phen, colour =spp_phen)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_continuous(limits = c(120,220), breaks = c(120,140,160,180,200,220)) +
  scale_x_continuous(limits = c(2001,2017), breaks = c(2005,2010,2015)) +
  theme_bw()

qhi_janet <- qhi_janet %>% mutate(spp_phen = paste0(spp, "_", phenophase))

SALARC_means_janet_plot <- ggplot(qhi_janet[qhi_janet$spp == "SALARC",], aes(x = year, y = day, group= spp_phen, colour =spp_phen)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_continuous(limits = c(140,250), breaks = c(140,160,180,200,220,240)) +
  scale_x_continuous(limits = c(2001,2017), breaks = c(2005,2010,2015)) +
  theme_bw()
ERIVAG_means_janet_plot <- ggplot(qhi_janet[qhi_janet$spp == "ERIVAG",], aes(x = year, y = day, group= spp_phen, colour =spp_phen)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_continuous(limits = c(130,170), breaks = c(130,140,150,160, 170)) +
  scale_x_continuous(limits = c(2001,2017), breaks = c(2005,2010,2015)) +
  theme_bw()
DRYINT_means_janet_plot <- ggplot(qhi_janet[qhi_janet$spp == "DRYINT",], aes(x = year, y = day, group= spp_phen, colour =spp_phen)) +
  geom_point() +
  geom_smooth(method = lm) +
  scale_y_continuous(limits = c(120,220), breaks = c(120,140,160,180,200,220)) +
  scale_x_continuous(limits = c(2001,2017), breaks = c(2005,2010,2015)) +
  theme_bw()

library(gridExtra)
panel_layout <- rbind(c(1,2,3),c(4,5,6))
png(filename = paste0(script_path, "/dataset_compar.png"), width = 1000, height = 600)
grid.arrange(SALARC_means_plot,
             ERIVAG_means_plot,
             DRYINT_means_plot,
             SALARC_means_janet_plot,
             ERIVAG_means_janet_plot,
             DRYINT_means_janet_plot,
             layout_matrix = panel_layout)
dev.off()