# Quick quality control script to assess effect of using the wrong cell sizes
# Jakob Assmann j.assmann@bio.au.dk 25 Jan 2021

# Dependencies
library(tidyverse)
library(cowplot)

# Load data
load("data/coastal_phen_new.Rda")
coastal_phen_new <- coastal_phen
load("data/coastal_phen.Rda")

# Check structure of the data frames is the same
ncol(coastal_phen) == ncol(coastal_phen_new)
nrow(coastal_phen) == nrow(coastal_phen_new)
names(coastal_phen) == names(coastal_phen_new)

# Compare difference by column
lapply(1:ncol(coastal_phen), function(x){
  return(c(names(coastal_phen)[x], identical(coastal_phen[,x], coastal_phen_new[,x])))
})

# Okay bad news, ice colums are all different
# calculate difference
onset_ice_diff <- coastal_phen$onset_ice_melt - coastal_phen_new$onset_ice_melt
hist(onset_ice_diff)

# Look at what this means on a site and year basis
ice_site_year <- coastal_phen %>% distinct(site_name, year, onset_ice_melt)
ice_site_year_new <- coastal_phen_new %>% distinct(site_name, year, onset_ice_melt)

# Plot differences
diff_plot <- ggplot() +
  geom_line(aes(x = year, y = onset_ice_melt), ice_site_year_new, col = "red") +
  geom_point(aes(x = year, y = onset_ice_melt), ice_site_year_new, col = "red") +
  geom_ribbon(aes(x = ice_site_year$year, 
              ymin = ice_site_year$onset_ice_melt,
              ymax = ice_site_year_new$onset_ice_melt), 
              ice_site_year,
              fill = "red") +
  geom_line(aes(x = year, y = onset_ice_melt), ice_site_year, col = "blue") +
  geom_point(aes(x = year, y = onset_ice_melt), ice_site_year, col = "blue") +
  facet_wrap(vars(site_name), scales = "free", ) +
  labs(x = "Year",  y = "Spring drop sea-ice (DoY)") +
  theme_cowplot() +
  theme(strip.background = element_rect(fill = "white"))
save_plot("data/2020_02_01_quality_control/qc_2021/diff_sea_ice_drop.png",
          diff_plot,
          base_height = 8)
