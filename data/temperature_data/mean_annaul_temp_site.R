# Mean Annual temperature change plots for the study sites in the 
# coastal phenology dataset.

# load dependencies
library(dplyr)
library(ggplot2)
library(cowplot)

# load datasets
alexfiord_temp <- read.csv("scripts/users/jassmann/phenology/temperature_data/alexfiord/alexfiord_daily_temp.csv")
barrow_temp <- read.csv("scripts/users/jassmann/phenology/temperature_data/barrow/barrow_daily_temp.csv")
qhi_temp <- read.csv("scripts/users/jassmann/phenology/temperature_data/qhi/qhi_daily_temp.csv")
zackenberg_temp <- read.csv("scripts/users/jassmann/phenology/temperature_data/zackenberg/zackenberg_daily_temp.csv")

temp_data <- bind_rows(alexfiord_temp, 
                       barrow_temp, 
                       qhi_temp, 
                       zackenberg_temp)  %>% group_by(site_name)

temp_data %>% summarise(min_year = min(year),
                        max_year = max(year))

temp_data <- temp_data %>% 
  filter(year >= 1995, year <= 2015) %>% 
  group_by(site_name, year) %>% summarise(mat = mean(temp, na.rm = T))

temp_trends <- lapply(unique(temp_data$site_name), function(site){
  site_temp <- filter(temp_data, site_name == site)
  site_lm  <- lm(mat ~ year, site_temp)
  print(site)
  print(summary(site_lm))
  site_summary <- summary(site_lm)
  site_slope <- site_summary$coefficients[[2]]
  return(setNames(c(site, site_slope), c("site_name", "slope")))
})
temp_trends <- do.call(rbind.data.frame, temp_trends)
names(temp_trends) <- c("site_name", "slope")
temp_trends$slope <- as.numeric(as.character(temp_trends$slope))

ylab_temp <- "Mean Annual Temperature (°C)"

site_names <- unique(as.character(temp_data$site_name))
colour_theme_sites <- data.frame(site_name = c(site_names),
                                 name_pretty =  c("Alexandra Fiord", 
                                                  "Utqiaġvik",
                                                  "Qikiqtaruk",
                                                  "Zackenberg"),
                                 colour = c("#324D5CFF", 
                                            "#46B29DFF", 
                                            "#C2A33EFF",
                                            "#E37B40FF"),
                                 stringsAsFactors = F) # 5th colour #F53855FF
factor(temp_data$site_name)
temp_plot <- ggplot(temp_data, aes(x = year, y = mat, group = site_name, colour = site_name)) + 
  geom_smooth(method = "lm") +
  ylab(ylab_temp) +
  xlab("") +
  scale_colour_manual(values = colour_theme_sites$colour) +
  annotate("text", x = 1995, y = -6, 
           label = paste0(colour_theme_sites$name_pretty[1], 
                          ": "),
           hjust = 0, colour = colour_theme_sites$colour[1])  +
  annotate("text", x = 1995, y = -6.5, 
           label = paste0(colour_theme_sites$name_pretty[2], 
                          ": "),
           hjust = 0, colour = colour_theme_sites$colour[2])  +
  annotate("text", x = 1995, y = -7, 
           label = paste0(colour_theme_sites$name_pretty[3], 
                          ": "),
           hjust = 0, colour = colour_theme_sites$colour[3])  +
  annotate("text", x = 1995, y = -7.5, 
           label = paste0(colour_theme_sites$name_pretty[4], 
                          ": "),
           hjust = 0, colour = colour_theme_sites$colour[4])  +
  
  annotate("text", x = 2002, y = -6, 
           label = paste0("+", round(temp_trends$slope[1] * 10, 1), "°C / decade"),
           hjust = 0, colour = colour_theme_sites$colour[1])  +
  annotate("text", x = 2002, y = -6.5, 
           label = paste0("+", round(temp_trends$slope[2] * 10, 1), "°C / decade*"),
           hjust = 0, colour = colour_theme_sites$colour[2])  +
  annotate("text", x = 2002, y = -7, 
           label = paste0("+", round(temp_trends$slope[3] * 10, 1), "°C / decade*"),
           hjust = 0, colour = colour_theme_sites$colour[3])  +
  annotate("text", x = 2002, y = -7.5, 
           label = paste0("+", round(temp_trends$slope[4] * 10, 1), "°C / decade*"),
           hjust = 0, colour = colour_theme_sites$colour[4])  +
  theme(legend.position = "none") 
save_plot("scripts/users/jassmann/phenology/temperature_data/mat_plot.png", temp_plot,
          base_aspect_ratio = 1.3 )

