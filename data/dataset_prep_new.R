### Coastal phenology prep script
# Jakob Assmann, j.assmann@ed.ac.uk February 2018

# Dependencies
library(dplyr)
library(ggplot2)

# Housekeeping
script_path <- "data/"
sites <- c("alexfiord", "barrow", "qhi", "zackenberg") # NB lower case for objects

# Import phenology data
alexfiord_phen <- read.csv("data/phenology_data/alexfiord/alexfiord_phen.csv")
barrow_phen <- read.csv("data/phenology_data/barrow/barrow_phen.csv")
qhi_phen <- read.csv("data/phenology_data/qhi/qhi_phen.csv")
zackenberg_phen <- read.csv("data/phenology_data/zackenberg/zackenberg_phen.csv")

# Import temperature data
alexfiord_temp <- read.csv("data/temperature_data/alexfiord/alexfiord_daily_temp.csv")
barrow_temp <- read.csv("data/temperature_data/barrow/barrow_daily_temp.csv")
qhi_temp <- read.csv("data/temperature_data/qhi/qhi_daily_temp.csv")
zackenberg_temp <- read.csv("data/temperature_data/zackenberg/zackenberg_daily_temp.csv")

# Import sea ice data
alexfiord_sea_ice <- read.csv("data/sea_ice_data/alexfiord/alexfiord_sea_ice_extent_new.csv")
barrow_sea_ice <- read.csv("data/sea_ice_data/barrow/barrow_sea_ice_extent_new.csv")
qhi_sea_ice <- read.csv("data/sea_ice_data/qhi/qhi_sea_ice_new.csv")
zackenberg_sea_ice <- read.csv("data/sea_ice_data/zackenberg/zackenberg_sea_ice_extent_new.csv")

### Prepare temperature data
extract_temperature <- function(site_name) {
  ## Housekeeping
  # create temporary dataframes for handling
  phen_data <- get(paste0(site_name, "_phen"))
  temp_data <- get(paste0(site_name, "_temp"))
  # get list of species, phenology stage combos
  spp_phen_combos <- unique(phen_data$spp_phen)
  # prep global dataframe for outputs
  phen_data$phase_temp <- NA
  phen_data$phase_temp_unmodified <- NA
  phen_data$site_phase_temp <- NA
  phen_data$site_phase_temp_unmodified <- NA
  phen_data$spring_temp <- NA
  phen_data$apr_temp <- NA
  phen_data$may_temp <- NA
  phen_data$jun_temp <- NA
  
  
  ## Extract time period dates for each spp_phen combo
  # initiate storage df
  time_periods_df <- data.frame(spp_phen = spp_phen_combos, 
                                start_doy = NA,
                                start_doy_unmodified = NA,
                                end_doy = NA
                                )
  # determine time period
  sapply(spp_phen_combos, function(spp_phen_comb){
    # subset data for easier working
    phen_data_subset <- phen_data %>% filter(spp_phen == spp_phen_comb)
    
    # determine start of period by extracting first day of snowmelt 
    # on record and subtracting 14 days (2 weeks)
    start_of_period <- min(phen_data_subset$snowmelt, na.rm = T) - 14
    # and now without the 14 day modifier
    start_of_period_unmodified <- min(phen_data_subset$snowmelt, na.rm = T)

    # calculate doy when 75% of phenology events
    # first extract leafout dates and sort ascending
    phen_doy <- phen_data_subset$doy %>% sort()
    # extract doy
    end_of_period <- phen_doy[round(length(phen_doy)*0.75)]
    
    # store in output df
    time_periods_df[time_periods_df$spp_phen == spp_phen_comb,]$start_doy <<- start_of_period
    time_periods_df[time_periods_df$spp_phen == spp_phen_comb,]$start_doy_unmodified <<- start_of_period_unmodified
    time_periods_df[time_periods_df$spp_phen == spp_phen_comb,]$end_doy <<- end_of_period
  })

  # Determine average daily temp from 14 days prior snowmelt to 75% of phenology at site level
  # 2 weeks prior snowmelt
  start_period_site <- min(phen_data$snowmelt, na.rm = T) - 14
  # wihtout modifier
  start_period_site_unmodified <- min(phen_data$snowmelt, na.rm = T)
  # Extract all observation doys on record and sort ascending
  site_phen_doy <- phen_data$doy %>% sort()
  # Determine 75% occurence doy
  end_period_site <- site_phen_doy[round(length(site_phen_doy)*0.75)]
  # Calculate site mean daily temp for the periods
  site_phase_temp <- temp_data %>% group_by(year) %>% 
    filter(doy >= start_period_site & doy <= end_period_site) %>%
    summarise(site_phase_temp = round(mean(temp, na.rm = T), 1))
  site_phase_temp_unmodified <- temp_data %>% group_by(year) %>% 
    filter(doy >= start_period_site_unmodified & doy <= end_period_site) %>%
    summarise(site_phase_temp_unmodified = round(mean(temp, na.rm = T), 1))
 
  # Save time periods data frame for QC purposes.
  write.csv(cbind(time_periods_df, start_period_site, start_period_site_unmodified, end_period_site), file= paste0(script_path, "quality_control/time_periods_temperature_", site_name, ".csv"))

  ## Calculate yearly averages for each time period and export to phenology dataframe
  sapply(time_periods_df$spp_phen, function(spp_phen_comb){
    # grab years on record and start / end dates
    years_on_record <- unique(phen_data[phen_data$spp_phen == spp_phen_comb,]$year)
    start_doy <- time_periods_df[time_periods_df$spp_phen == spp_phen_comb,]$start_doy
    start_doy_unmodified <- time_periods_df[time_periods_df$spp_phen == spp_phen_comb,]$start_doy_unmodified
    end_doy <- time_periods_df[time_periods_df$spp_phen == spp_phen_comb,]$end_doy
    # Calculate average temperatures for time period
    phase_temp <- temp_data %>% group_by(year) %>% 
      filter(doy >= start_doy & doy <= end_doy) %>%
      summarise(phase_temp = round(mean(temp, na.rm = T), 1))
    # and without 14 day modifier
    phase_temp_unmodified <- temp_data %>% group_by(year) %>% 
      filter(doy >= start_doy_unmodified & doy <= end_doy) %>%
      summarise(phase_temp_unmodified = round(mean(temp, na.rm = T), 1))

    # Calculate monthly and spring (april - june) averages
    monthly_temp <- temp_data %>% group_by(year, month) %>%
      summarise(temp = round(mean(temp, na.rm = T),1))
    spring_temp <- temp_data %>% group_by(year) %>%
      filter(month >= 5 & month <= 7) %>%
      summarise(temp = round(mean(temp, na.rm = T),1))

    # Export to phen_data df
    sapply(years_on_record, function(year_to_set){
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$phase_temp <<- phase_temp[phase_temp$year == year_to_set,]$phase_temp
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$phase_temp_unmodified <<- phase_temp_unmodified[phase_temp_unmodified$year == year_to_set,]$phase_temp_unmodified
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$site_phase_temp <<- site_phase_temp[site_phase_temp$year == year_to_set,]$site_phase_temp
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$site_phase_temp_unmodified <<- site_phase_temp_unmodified[site_phase_temp_unmodified$year == year_to_set,]$site_phase_temp_unmodified
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$spring_temp <<- spring_temp[spring_temp$year == year_to_set,]$temp
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$apr_temp <<- monthly_temp[monthly_temp$year == year_to_set & monthly_temp$month == 4,]$temp
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$may_temp <<- monthly_temp[monthly_temp$year == year_to_set & monthly_temp$month == 5,]$temp
      phen_data[phen_data$year == year_to_set & phen_data$spp_phen == spp_phen_comb,]$jun_temp <<- monthly_temp[monthly_temp$year == year_to_set & monthly_temp$month == 6,]$temp
    }) 
  })
  return(phen_data)
}

# Extract temperatures for each site and append to _phen dataframe
alexfiord_phen <- extract_temperature(sites[1])
barrow_phen <- extract_temperature(sites[2])
qhi_phen <- extract_temperature(sites[3])
zackenberg_phen <- extract_temperature(sites[4])

#### Next up calculate onset of sea ice melt and average spring extent
extract_sea_ice <- function(site_name){
  # House keeping
  phen_data <- get(paste0(site_name, "_phen"))
  sea_ice_data <- get(paste0(site_name, "_sea_ice"))
  sea_ice_data$date <- as.Date(as.character(sea_ice_data$date))
  
  # prep output collumns
  phen_data$onset_ice_melt <- NA
  phen_data$spring_extent <- NA
  phen_data$phase_extent <- NA
  phen_data$phase_extent_unmodified <- NA
  
  # determine max and 85% of sea ice extent
  max_extent <- max(sea_ice_data$sea_ice_extent, na.rm = T)
  extent_85 <- max_extent * 0.85
  
  # determine earliest date when annual minimum is reached
  onset_sea_ice_melt <- sea_ice_data %>% group_by(year) %>% summarise(min_extent = min(sea_ice_extent, na.rm = T))
  onset_sea_ice_melt$min_doy <- sapply(onset_sea_ice_melt$year, function(year){
    min_doy <- min(sea_ice_data[sea_ice_data$year == year &
                                  sea_ice_data$sea_ice_extent == onset_sea_ice_melt[
                                                                  onset_sea_ice_melt$year == year,]$min_extent
                                ,]$doy, na.rm = T)
  })
  
  # find last day at which the 85 extend is esceeded prior maximum is reached and add 1
  onset_sea_ice_melt$onset_melt <- sapply(onset_sea_ice_melt$year, function(year){
    onset_melt <- max(sea_ice_data[sea_ice_data$year == year &
                                     sea_ice_data$doy < onset_sea_ice_melt[onset_sea_ice_melt$year == year,]$min_doy &
                                     sea_ice_data$sea_ice_extent > extent_85,]$doy, na.rm = T) + 1
  })
    
  onset_sea_ice_melt <- onset_sea_ice_melt %>% mutate(date = as.Date(paste0(year, "-", onset_melt), format = "%Y - %j"))
  onset_sea_ice_melt$ice_extent <- sea_ice_data[sea_ice_data$date %in% onset_sea_ice_melt$date,]$sea_ice_extent
  # avoid NAs by setting missing to extent_85
  if(sum(is.na(onset_sea_ice_melt$ice_extent)) != 0){ 
    onset_sea_ice_melt[is.na(onset_sea_ice_melt$ice_extent),]$ice_extent <- round(extent_85)
  }
  
  if(site_name == "alexfiord"){
    site_name_pretty <- "Alexandra Fiord"
  } else if (site_name == "barrow") {
    site_name_pretty <- "Utqiaġvik"
  } else if (site_name == "qhi"){
    site_name_pretty <- "Qikiqtaruk"
  } else if(site_name == "zackenberg"){
    site_name_pretty <- "Zackenberg"
  } else {
    site_name_pretty <- "something went wrong"
  }
  # plot for quality control on sea ice trajectory
  ggplot(sea_ice_data, aes(x = date, y = sea_ice_extent)) +
    geom_line() + 
    ggtitle(label = paste0(site_name_pretty)) +
    scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
    geom_point(data = onset_sea_ice_melt, aes(x = date, y = ice_extent), colour = "firebrick", size = 3) +
    scale_y_continuous(limits = c(0, 150000), breaks = seq(0, 150000, 25000)) +
    ylab(expression(paste("Regional Sea Ice Extent (", km^2, ")"))) +
    xlab("") +
    theme_bw() + theme(legend.position = "none", 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.text.x = element_text(angle = 90, vjust = 1.7), 
                       axis.text = element_text(size = 20),
                       panel.border = element_rect(colour = "black"),
                       axis.title = element_text(size = 30),
                       plot.title = element_text(size = 30)) 
  ggsave(paste0(script_path, "sea_ice_data/", site_name, "/", site_name, "_sea_ice_new.png"), width = 16, height = 8)
 
  # spring (may june july) average
  onset_sea_ice_melt$spring_extent <- sea_ice_data %>%
    group_by(year) %>% filter(month >= 5 & month <= 7) %>%
    summarise(spring_avg = round(mean(sea_ice_extent, na.rm = T))) %>% select(spring_avg) %>% t() %>% as.vector()
  # convert into percentage with 1 decimal
  onset_sea_ice_melt$spring_extent  <- round(onset_sea_ice_melt$spring_extent / max_extent * 100, 1)
  
  # Determine phen phase dates for the site
  # 2 weeks prior snowmelt
  start_period_site <- min(phen_data$snowmelt, na.rm = T) - 14
  # wihtout modifier
  start_period_site_unmodified <- min(phen_data$snowmelt, na.rm = T)
  # Extract all observation doys on record and sort ascending
  site_phen_doy <- phen_data$doy %>% sort()
  # Determine 75% occurence doy
  end_period_site <- site_phen_doy[round(length(site_phen_doy)*0.75)]
  # write dates for QC purposes
  write.csv(cbind(start_period_site, start_period_site_unmodified, end_period_site), file = paste0(script_path, "/quality_control/time_period_ice_", site_name, ".csv"))
  
  # spring phen phase average
  onset_sea_ice_melt$phase_extent <- sea_ice_data %>%
    group_by(year) %>% filter(doy >= start_period_site & doy <= end_period_site) %>%
    summarise(phase_avg = round(mean(sea_ice_extent, na.rm = T))) %>% select(phase_avg) %>% t() %>% as.vector()
  # convert into percentage with 1 decimal
  onset_sea_ice_melt$phase_extent  <- round(onset_sea_ice_melt$phase_extent / max_extent * 100, 1)
  # spring phen phase average unmodified
  onset_sea_ice_melt$phase_extent_unmodified <- sea_ice_data %>%
    group_by(year) %>% filter(doy >= start_period_site_unmodified & doy <= end_period_site) %>%
    summarise(phase_unmodified_avg = round(mean(sea_ice_extent, na.rm = T))) %>% select(phase_unmodified_avg) %>% t() %>% as.vector()
  # convert into percentage with 1 decimal
  onset_sea_ice_melt$phase_extent_unmodified  <- round(onset_sea_ice_melt$phase_extent_unmodified / max_extent * 100, 1)

  
  # export to phen_data
  years_on_record <- unique(phen_data$year)
  sapply(years_on_record, function(year_to_set){
    phen_data[phen_data$year == year_to_set,]$onset_ice_melt <<- onset_sea_ice_melt[onset_sea_ice_melt$year == year_to_set,]$onset_melt
   
    phen_data[phen_data$year == year_to_set,]$spring_extent <<- onset_sea_ice_melt[onset_sea_ice_melt$year == year_to_set,]$spring_extent
    phen_data[phen_data$year == year_to_set,]$phase_extent <<- onset_sea_ice_melt[onset_sea_ice_melt$year == year_to_set,]$phase_extent
    phen_data[phen_data$year == year_to_set,]$phase_extent_unmodified <<- onset_sea_ice_melt[onset_sea_ice_melt$year == year_to_set,]$phase_extent_unmodified
    
  }) 
  
  # return dataframe
  return(phen_data)
}


# no sea_ice data avialable for 2017 throw out 2017 data from sea_ice dataset
# that is not needed
alexfiord_sea_ice <- alexfiord_sea_ice[alexfiord_sea_ice$year >= 1992 & alexfiord_sea_ice$year <= 2013,]
barrow_sea_ice <- barrow_sea_ice[barrow_sea_ice$year >= 1994 & barrow_sea_ice$year <= 2014,]
qhi_sea_ice <- qhi_sea_ice[qhi_sea_ice$year >= 2001 & qhi_sea_ice$year <= 2016,]
zackenberg_sea_ice <- zackenberg_sea_ice[zackenberg_sea_ice$year >= 1996 & zackenberg_sea_ice$year <= 2011,]

# no sea_ice data available for 2017 throw out QHI data
qhi_phen <- qhi_phen[qhi_phen$year != 2017,]

# map extract_sea_ice function
alexfiord_phen <- extract_sea_ice(sites[1])
barrow_phen <-  extract_sea_ice(sites[2])
qhi_phen <-  extract_sea_ice(sites[3])
zackenberg_phen <-  extract_sea_ice(sites[4])

## Merge data frames
coastal_phen <- rbind(alexfiord_phen,
                      barrow_phen,
                      qhi_phen,
                      zackenberg_phen)

## create year as factor collumn
coastal_phen$year_fac <- factor(coastal_phen$year)
# re arrange
coastal_phen <- coastal_phen[c(1,2,3,4,5,6,7,23,seq(8,22))]

# export data
save(coastal_phen, file = paste0(script_path, "coastal_phen_new.Rda"))
write.csv(coastal_phen, file = paste0(script_path, "coastal_phen_new.csv"), row.names = F)

# Plot mean of all variables for each phenostage for QC
mean_phenology <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_phenology = mean(doy, na.rm = T))
ggplot(data = mean_phenology, aes(x = site_spp_phen, y = mean_phenology)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_phenology_new.png"), width = 16, height = 8)

mean_snowmelt <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_snowmelt = mean(snowmelt, na.rm = T))
ggplot(data = mean_snowmelt, aes(x = site_spp_phen, y = mean_snowmelt)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_snowmelt_new.png"), width = 16, height = 8)

mean_phase_temp <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_phase_temp = mean(phase_temp, na.rm = T))
ggplot(data = mean_phase_temp, aes(x = site_spp_phen, y = mean_phase_temp)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_phase_temp_new.png"), width = 16, height = 8)

mean_phase_temp_unmodified <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_phase_temp_unmodified = mean(phase_temp_unmodified, na.rm = T))
ggplot(data = mean_phase_temp_unmodified, aes(x = site_spp_phen, y = mean_phase_temp_unmodified)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_phase_temp_unmodified_new.png"), width = 16, height = 8)

mean_site_phase_temp <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_site_phase_temp = mean(phase_temp, na.rm = T))
ggplot(data = mean_site_phase_temp, aes(x = site_spp_phen, y = mean_site_phase_temp)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_site_phase_temp_new.png"), width = 16, height = 8)

mean_site_phase_temp_unmodified <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_site_phase_temp_unmodified = mean(phase_temp_unmodified, na.rm = T))
ggplot(data = mean_site_phase_temp_unmodified, aes(x = site_spp_phen, y = mean_site_phase_temp_unmodified)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_site_phase_temp_unmodified_new.png"), width = 16, height = 8)

mean_onset_ice_melt <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_onset_ice_melt = mean(onset_ice_melt, na.rm = T))
ggplot(data = mean_onset_ice_melt, aes(x = site_spp_phen, y = mean_onset_ice_melt)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_onset_ice_melt_new.png"), width = 16, height = 8)

mean_spring_extent <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_spring_extent = mean(spring_extent, na.rm = T))
ggplot(data = mean_spring_extent, aes(x = site_spp_phen, y = mean_spring_extent)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_spring_extent_new.png"), width = 16, height = 8)

mean_phase_extent <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_phase_extent = mean(phase_extent, na.rm = T))
ggplot(data = mean_phase_extent, aes(x = site_spp_phen, y = mean_phase_extent)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_phase_extent_new.png"), width = 16, height = 8)

mean_phase_extent_unmodified <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(mean_phase_extent_unmodified = mean(phase_extent_unmodified, na.rm = T))
ggplot(data = mean_phase_extent_unmodified, aes(x = site_spp_phen, y = mean_phase_extent_unmodified)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/mean_phase_extent_unmodified_new.png"), width = 16, height = 8)

phenology <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = phenology, aes(x = site_spp_phen, y = doy)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/phenology_new.png"), width = 16, height = 8)

snowmelt <- coastal_phen %>% group_by(site_spp_phen)
ggplot(data = snowmelt, aes(x = site_spp_phen, y = snowmelt)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/snowmelt_new.png"), width = 16, height = 8)

phase_temp <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = phase_temp, aes(x = site_spp_phen, y = phase_temp)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/phase_temp_new.png"), width = 16, height = 8)

phase_temp <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = phase_temp, aes(x = site_spp_phen, y = spring_temp)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/spring_temp_new.png"), width = 16, height = 8)

phase_temp_unmodified <- coastal_phen %>% group_by(site_spp_phen)
ggplot(data = phase_temp_unmodified, aes(x = site_spp_phen, y = phase_temp_unmodified)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/phase_temp_unmodified_new.png"), width = 16, height = 8)

site_phase_temp <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = site_phase_temp, aes(x = site_spp_phen, y = site_phase_temp)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/site_phase_temp_new.png"), width = 16, height = 8)

site_phase_temp_unmodified <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = site_phase_temp_unmodified, aes(x = site_spp_phen, y = site_phase_temp_unmodified)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/site_phase_temp_unmodified_new.png"), width = 16, height = 8)

onset_ice_melt <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = onset_ice_melt, aes(x = site_spp_phen, y = onset_ice_melt)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/onset_ice_melt_new.png"), width = 16, height = 8)

spring_extent <- coastal_phen %>% group_by(site_spp_phen)
ggplot(data = spring_extent, aes(x = site_spp_phen, y = spring_extent)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/spring_extent_new.png"), width = 16, height = 8)

phase_extent <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = phase_extent, aes(x = site_spp_phen, y = phase_extent)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/phase_extent_new.png"), width = 16, height = 8)

phase_extent_unmodified <- coastal_phen %>% group_by(site_spp_phen) 
ggplot(data = phase_extent_unmodified, aes(x = site_spp_phen, y = phase_extent_unmodified)) + geom_point() + theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(script_path, "quality_control/phase_extent_unmodified_new.png"), width = 16, height = 8)


