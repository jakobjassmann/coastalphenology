## Script to generate daily means of alexfiord temperature data

# dependencies
library(dplyr)

# set script path
script_path <- "scripts/users/jassmann/phenology/temperature_data/alexfiord/"

# load temperature file
alexfiord_temp <- read.csv("scripts/users/jassmann/phenology/temperature_data/alexfiord/AirTforAnne_Thru2015_gapfilled.csv")[-1]
names(alexfiord_temp) <- c("year", "doy", "temp")

alexfiord_temp <- alexfiord_temp %>% 
  mutate(date = as.Date(paste0(year, "-", doy), format = "%Y-%j")) %>%
  mutate(month = format(date, "%m"), day = format(date, "%d"))
alexfiord_temp$site_name <- "ALEXFIORD"
names(alexfiord_temp)
alexfiord_temp <- alexfiord_temp[,c(7,4,1,5,6,2,3)]

# quality check!
ggplot(alexfiord_temp, aes(x = date, y = temp)) + geom_point() + theme_bw()


# save
save(alexfiord_temp, file = paste0(script_path, "alexfiord_daily_temp.Rda"))
write.csv(alexfiord_temp, file = paste0(script_path, "alexfiord_daily_temp.csv"), row.names=FALSE)
