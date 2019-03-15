## Script to generate daily means of barrow temperature data

# dependencies
library(dplyr)

# set script path
script_path <- "scripts/users/jassmann/phenology/temperature_data/zackenberg/"

# load temperature file
zackenberg_temp <- read.csv("scripts/users/jassmann/phenology/temperature_data/zackenberg/View_ClimateBasis_Zackenberg_Data_Temperature_Air_temperature_200_cm__60min_average_deg_C120220181756020442.csv")

# tidyp up names
names(zackenberg_temp) <- c("date", "time", "hourly_temp") 

# turn date column into date
zackenberg_temp$date = as.Date(zackenberg_temp$date)

# chuck out absent readings
zackenberg_temp$hourly_temp[zackenberg_temp$hourly_temp == -9999.0]<- NA

zackenberg_daily_temp <- zackenberg_temp %>% 
  mutate(year = format(date, "%Y"), month = format(date, "%m"), day = format(date, "%d"), doy = format(date, "%j")) %>%
  group_by(year, month, day) %>% 
  summarise(temp = round(mean(hourly_temp), 1)) %>% 
  mutate(date = as.Date(paste0(year, "-", month, "-", day))) %>%
  mutate(doy = format(date, "%j"))
zackenberg_daily_temp$site_name <- "ZACKENBERG"
names(zackenberg_daily_temp)
zackenberg_daily_temp <- zackenberg_daily_temp[c(7,5,1,2,3,6,4)]
names(zackenberg_daily_temp)

# quality check!
ggplot(zackenberg_daily_temp, aes(x = date, y = temp)) + geom_point() + theme_bw()

# save
save(zackenberg_daily_temp, file = paste0(script_path, "zackenberg_daily_temp.Rda"))
write.csv(zackenberg_daily_temp, file = paste0(script_path, "zackenberg_daily_temp.csv"), row.names = F)
