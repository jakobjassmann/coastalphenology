### Script to extract the Environment Canada daily weather station data
### Jakob Assmann 8 November 2016

require(dplyr)

# First lets create a vector of all the weather station files contained in the folder of interest
# for this we use a variable containg the folder path of interest and the Sys.glob command for
# globbing the file paths (wildcard expansion on file paths)
folder_path <- "scripts/users/jassmann/phenology/temperature_data/qhi/"
files <- Sys.glob(paste0(folder_path, "eng*"))

# check that it has worked:
files
# nice!

# read in all files and bind them into one data frame bind_rows() is super helpful!
qhi_weather <- files %>% lapply(read.csv, skip = 25, stringsAsFactors = FALSE) %>% lapply(select, Date.Time, Year, Month, Day, Data.Quality, Mean.Temp...C., Mean.Temp.Flag) %>% bind_rows()

# excellent now we have one big data frame of all the years. time to get familiar with it.

# collumn names
names(qhi_weather)

# actually we are only interested in the mean temeprature so lets get rid of the other stuff
qhi_weather <- qhi_weather %>% select(Date.Time, Year, Month, Day, Data.Quality, Mean.Temp...C., Mean.Temp.Flag)

# Get rows that have a temperature flag
flagged_rows <- qhi_weather[which(qhi_weather$Mean.Temp.Flag != ""),]
# let's check uot what the flags are
unique(flagged_rows$Mean.Temp.Flag)
# It turns out it's just missing (M) and estimated values (E)

# Get rows with estimated values
estimated <- which(qhi_weather$Mean.Temp.Flag == "E")
# Damn, there are 83 estimated values. Actually that is not that much. Let's just leave them there
# qhi_weather[estimated,6] <- NA

# Get rows with missing values 
missing <- which(qhi_weather$Mean.Temp.Flag == "M")

# Let's chuck out all estimated and 'missing' values (weirdly some are marked missing but have a value - strange!)
qhi_weather[which(qhi_weather$Mean.Temp.Flag != ""),]$Mean.Temp...C. <- NA

# Get rid of all but the date and mean temp columns
qhi_weather <- qhi_weather %>% select(-Data.Quality, -Mean.Temp.Flag)

# Now, let's find out how many values are missing
nrow(qhi_weather[is.na(qhi_weather$Mean.Temp...C.),])
# of total
nrow(qhi_weather)
# 1775 of 5844 = 30% quite a lot. we will have to do some gap filling.
# maybe another day...

# let's clean up first:
rm(flagged_rows, estimated, missing)

### the next part of the script requires the ../komakuk/extract_weather_data script to be run!

# let's combine the two into one dataset
qhi_phen_temp <- data.frame(qhi_weather[,1:4], qhi_temp = qhi_weather$Mean.Temp...C., kom_temp =komakuk_weather$Mean.Temp...C.)

# and now some tidying up and convert dates into year + doy
# first turn dates from characters into formated dates
qhi_phen_temp$Date.Time <- as.Date(qhi_phen_temp$Date.Time)
# now tidy up and creat doy collum
qhi_phen_temp <- qhi_phen_temp %>% mutate(doy = as.numeric(format(Date.Time, "%j"))) %>% select(Year, Month, Day, doy, qhi_temp, kom_temp)

# let's see how they correlate
cor(qhi_phen_temp$qhi_temp, qhi_phen_temp$kom_temp, use = "na.or.complete")
# 0.9864276 that is a pretty good correlation.

# create let's model qhi date based on kom data
qhi_lm <- lm(qhi_phen_temp$qhi_temp ~ qhi_phen_temp$kom_temp) 
summary(qhi_lm)
# looks good adjusted R^2 is 0.973 very high explanatory capabillity
# extract coefficients
qhi_lm_coef <- coef(qhi_lm)
# do a scatter plot with model for pretyness
plot(qhi_phen_temp$kom_temp, qhi_phen_temp$qhi_temp, xlab = "Komakuk Weather Station daily Temp (°C)", ylab = "QHI Weather Station daily Temp (°C)")
abline(qhi_lm_coef[1],qhi_lm_coef[2], col = "red")
text(10,-20, paste0("y = ", round(qhi_lm_coef[1], 3), " + ", round(qhi_lm_coef[2], 3)), col = "red")

# now use kom data to model qhi data
qhi_phen_temp <- qhi_phen_temp %>% mutate(qhi_mod_temp = (kom_temp*qhi_lm_coef[2] + qhi_lm_coef[1]))
# let's do a quick scatter plot to see how well the model compares to the actual values
# subset all values that are not NA
qhi_phen_temp_sub <- qhi_phen_temp %>% filter(!is.na(qhi_temp), !is.na(qhi_mod_temp))
plot(qhi_phen_temp_sub$qhi_mod_temp, qhi_phen_temp_sub$qhi_temp, xlab = "Modelled daily mean temp QHI (°C)", ylab = "Actual daily mean temp QHI (°C)")
abline(0,1, col = "red")
# looks good, tiday up
rm(qhi_phen_temp_sub)

# next complete the missing values of the qhi_temp with the modelled temperatrues
# best to do this by subsequent merging, first let's create a new collum
qhi_phen_temp <- data.frame(qhi_phen_temp, temp = qhi_phen_temp$qhi_temp)
# now we merge qhi_mod_temp in:
qhi_phen_temp$temp[is.na(qhi_phen_temp$qhi_temp)] <- qhi_phen_temp$qhi_mod_temp[is.na(qhi_phen_temp$qhi_temp)]
# this command is a bit of a nut to crack, but simply said it assigns values to all missing values in qhi_temp with those values in qhi_mod_temp that are missing in qhi_temp
# let's just round those temperature values to something meaningful
qhi_phen_temp$temp <- round(qhi_phen_temp$temp, 1)

# Quick stop in between. Let's check how many missing vlaues we still have.
sum(is.na(qhi_phen_temp$temp))
# 134 / 5844 that is not too bad! I'm happy to live with it

# do a quick monthly summary to check for monthly correlation (motivated by bad correlaiton with cru data of the final data)
qhi_phen_temp_mont <- qhi_phen_temp %>% group_by(Year, Month) %>% summarise(qhi_month = mean(qhi_temp), kom_month = mean(kom_temp), qhi_mod_month = mean(qhi_mod_temp))
cor(qhi_phen_temp_mont$qhi_month, qhi_phen_temp_mont$kom_month, use= "na.or.complete" )
# 0.997 => that is really not bad. How many months are there to correlate though?
sum(!is.na(qhi_phen_temp_mont$qhi_month) & !is.na(qhi_phen_temp_mont$kom_month))
# 33 not bad either
# do the same for spring and spring months (4-6)
qhi_phen_temp_summer_mont <- qhi_phen_temp %>% filter(Month >=4, Month <=6) %>% group_by(Year, Month) %>% summarise(qhi_month = mean(qhi_temp), kom_month = mean(kom_temp), qhi_mod_month = mean(qhi_mod_temp))
cor(qhi_phen_temp_summer_mont$qhi_month, qhi_phen_temp_summer_mont$kom_month, use= "na.or.complete" )
# 0.995 => that is really not bad. How many months are there to correlate though?
sum(!is.na(qhi_phen_temp_summer_mont$qhi_month) & !is.na(qhi_phen_temp_summer_mont$kom_month))
# 12 is not bad either


### Final prep and data export 

# tidy up
qhi_phen_temp <- qhi_phen_temp %>% select(-qhi_temp, -kom_temp, -qhi_mod_temp) %>%
  mutate(date = as.Date(paste0(Year, "-", Month, "-", Day)))
qhi_phen_temp$site_name <- "QHI"
qhi_daily_temp <- qhi_phen_temp[, c(7,6,1:5)]
names(qhi_daily_temp) <- tolower(names(qhi_daily_temp))

# Quality Check!

ggplot(qhi_daily_temp, aes(x = date, y = temp)) + geom_point() + theme_bw()

# Export Files
write.csv(qhi_daily_temp, "scripts/users/jassmann/phenology/temperature_data/qhi/qhi_daily_temp.csv", row.names = F)
save(qhi_daily_temp, file = "scripts/users/jassmann/phenology/temperature_data/qhi/qhi_daily_temp.Rda")

