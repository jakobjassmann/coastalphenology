## Script to generate daily means of barrow temperature data

# dependencies
library(dplyr)

# set script path
script_path <- "scripts/users/jassmann/phenology/temperature_data/barrow/"

# load temperature files

file_list <- list.files(path = script_path, pattern = "*.txt")

# create and fill adata frame with files
lapply(paste0(script_path, file_list), 
       function(x) {
         # assign(gsub(".*/met_brw_insitu_1_obop_hour_",
         #           "barrow_",
         #           gsub(".txt$", "", x)), 
         #      read.csv(x,sep = "", strip.white = T, col.names = c("site_code", "year", "month", "day", "hour", "wind_dir", "wind_speed", "wind_steadf", "pressure", "temp_2m", "temp_10m", "temp_top", "rel_hum", "precip")),
         #      envir = .GlobalEnv)
         
         # read in file
         temp_df <- read.csv(x,sep = "", strip.white = T, 
                             col.names = c("site_code", "year", "month", "day", 
                                           "hour", "wind_dir", "wind_speed", "wind_steadf", 
                                           "pressure", "temp_2m", "temp_10m", "temp_top", "rel_hum", "precip"))
         # check presence of barrow_all df
         if (exists("barrow_all")) {
           barrow_all <<- rbind(barrow_all, temp_df)
         } else {
           barrow_all <<- temp_df
         }
         return("Done")
       })

# calculate daily mean
barrow_daily_temp <- barrow_all %>% 
  select(year, month, day, hour, temp_2m) %>% 
  group_by(year, month, day) %>% 
  summarise(temp = round(mean(temp_2m), 1)) %>%
  mutate(date = as.Date(paste0(year, "-", month, "-", day))) %>%
  mutate(doy = format(date, "%j"))
barrow_daily_temp$site_name <- "BARROW"
barrow_daily_temp <- barrow_daily_temp[c(7,5,1,2,3,6,4)]

# quality check!
ggplot(barrow_daily_temp[barrow_daily_temp$month >= 4 & barrow_daily_temp$month <= 6,], aes(x = date, y = temp)) + geom_point() + theme_bw()
# quite a lot of erronous data, clean up! Chuck out anytihing that is below -45!
barrow_daily_temp[barrow_daily_temp$temp < -40,]$temp <- NA

# save
save(barrow_daily_temp, file = paste0(script_path, "barrow_daily_temp.Rda"))
write.csv(barrow_daily_temp, file = paste0(script_path, "barrow_daily_temp.csv"), row.names=FALSE )
