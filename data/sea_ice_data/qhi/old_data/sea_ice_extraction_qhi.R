### The purpose of this script is to extract Short Microwave Radiaiton Sea-Ice data from
# the NSCID Cilmate Recrod Data Set and prepare the data for inclusion in the phenology models
# Jakob Assmann (j.assmann@ed.ac.uk) August 2017

# load (some of the) dependencies
library(chron)
library(lattice)
library(ncdf4)
library(dplyr)
library(ggplot2)

################################################################################
### Extract Sea Ice Data from NSCID CDR

### Load sample Climate Data Record NetCDF file and familiarise
#####
# Set a file path to a NetCDF file
file_path <- "/Volumes/csce/biology/users/s1043792/sea_ice/sidads.colorado.edu/pub/DATASETS/NOAA/G02202_v2/north/daily/2001/seaice_conc_daily_nh_f13_20010501_v02r00.nc"

# load file (using ncdf4)
ncin <- nc_open(file_path)
print(ncin)

# let's look at longitude data
lon <- ncvar_get(ncin, "longitude")
nlon <- dim(lon)
head(lon)
lat <- ncvar_get(ncin, "latitude")
nlat <- dim(lat)
head(lat)
# They're in decimal degrees

## Lesson learned
# So these are big matrices with lat and long coordinates for each cell centre
# as the grid is a sterographic projection we do have individual lat/long for eachgrid cell and not just the margins as would only be required if the projection would not be warped
#####

### Identify the grid cell containing the field site and set bounding box (in this case QHI)
#####
# Set location
site_lat_lon <- c(69.57, -138.91) # Pauline Cove
# origin of the polar steographic grid is in the pole with positive x going twoards russia and negative y towards canada
# negative y is greenland and positive y is japan

# create a matrix containg the RMS difference between site_lat_long and grid_cell lat long
rms_distance <- matrix(nrow = 304, ncol = 448)
lat_diff <- 0
lon_diff <- 0
for(i in 1:304) {
  for(j in 1:448){
    lat_diff <- lat[i, j] - site_lat_lon[1]
    lon_diff <- lon[i, j] - site_lat_lon[2]
    rms_distance[i, j] <-  sqrt(0.5*(lat_diff^2 + lon_diff^2))
  }
}

# return cell array index of rms_distance minimum
centre_index <- which(rms_distance == min(rms_distance), arr.ind = T)
centre_index

###
# The cell contiaining QHI is: 65 and 228 !
###

# Confirm this by checking lat long ouptut
centre_lat_lon <- c(lat[centre_index], lon[centre_index])
centre_lat_lon

# Clauclate boundary box in xy grid that has the QHI cell in centre (going 10 rows/collumns in each direction)
x_range <- c(centre_index[1] - 10, centre_index[1] + 10)
y_range <- c(centre_index[2] - 10, centre_index[2] + 10)

#####

### Load pixel area data from binary file. 
#####
# prep matrix for storing data
pixel_area <- matrix(nrow = 304, ncol = 448)
# read in binary pixel area file
to_read <- file("/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/psn25area_v3.dat", "rb")
# read in pixel area values integer by integer
for (i in 1:304){
  for (j in 1:448){
    pixel_area[i, j] <- readBin(to_read, integer(), endian = "little", size = 4) / 1000 
    # NB the "endian" setting is really important, all other files provided by the NSCID with the pixel area file (masks etc) have a little endian, so we can assume that applies to the pixel area file too (though not specified). Also tested big endian to be sure and results do not make sense!
  }
}
# check whether there is anything left (if at end the value will be NULL)
readBin(to_read, integer(), endian = "little", size = 4)
# good!
# close file:
close(to_read)

# check whether all values got read:
pixel_area[304,448] # this is the final cell in the data set
# looks good!
#####

### Preparations for creating a mask
#####
# Calculate half way points along the grid for the selected cells within xrang and yrange
# these will allow us to recreate the grid lines in QGIS using the points2one polygon plug-in
# The grid lines can then be imported into QHIS/GoogleMaps/GoogleEarth to identify land and coastal pixels

# Prepare an empty data frame for the coordinates of cell corners of the x-grid (23 rows and 23 collumns, 21 for our bounding box and 2 for the outside lines)

grid_half_x <- data.frame(corner = rep(NA, (23*23)),
                          x_grid = rep(0, (23*23)),
                          y_grid = rep(0, (23*23)), 
                          latitude = rep(0, (23*23)),
                          longitude = rep(0, (23*23)))
# Fill the data frame

for(i in 1:23){
  for(j in 1:23){
    count <- (i*23-23) + (j-1) + 1 # this counter turns i and j into one continous count from 1 to 23*23 (529)
    x_corner <- ((x_range[1]-2) + i) # the left hand edge of the grid cell
    y_corner <- ((y_range[1]-2) + j) # the bottom edge of the grid cell
    grid_half_x$corner[count] <- paste0("x", x_corner, "y", y_corner)
    grid_half_x$x_grid[count] <- x_corner
    grid_half_x$y_grid[count] <- y_corner
    # now calculate lat long for each corner
    grid_half_x$latitude[count] <- 0.5* (lat[x_corner, y_corner] + lat[(x_corner + 1), (y_corner)]) # move along the x grid by 1
    grid_half_x$longitude[count] <- 0.5* (lon[x_corner, y_corner] + lon[(x_corner + 1), (y_corner)]) # move along the x grid by 1
  }
}


# Prepare an empty data frame for the coordinates of cell corners of the x-g

grid_half_y <- data.frame(corner = rep(NA, (23*23)), x_grid = rep(0, (23*23)), y_grid = rep(0, (23*23)), latitude = rep(0, (23*23)), longitude = rep(0, (23*23)))

for(i in 1:23){
  for(j in 1:23){
    count <- (i*23-23) + (j-1) +1
    x_corner <- ((x_range[1]-2) + i) # the left hand edge of the grid cell
    y_corner <- ((y_range[1]-2) + j) # the bottom edge of the grid cell
    grid_half_y$corner[count] <- paste0("x", x_corner, "y", y_corner)
    grid_half_y$x_grid[count] <- x_corner
    grid_half_y$y_grid[count] <- y_corner
    # now calculate lat long
    grid_half_y$latitude[count] <- 0.5* (lat[x_corner, y_corner] + lat[(x_corner), (y_corner + 1)]) # move along the y grid by 1
    grid_half_y$longitude[count] <- 0.5* (lon[x_corner, y_corner] + lon[(x_corner), (y_corner + 1)]) # move along the y grid by 1
  }
}

# export as CSV:
write.csv(grid_half_x,"/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/grid_half_x.csv")
write.csv(grid_half_y,"/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/grid_half_y.csv")

# also create file of grid centres
grid_centres <- data.frame(centre = rep(NA, (21*21)), x_grid = rep(0, (21*21)), y_grid = rep(0, (21*21)), latitude = rep(0, (21*21)), longitude = rep(0, (21*21)))
for(i in 1:21){
  for(j in 1:21){
    count <- (i*21-21) + (j-1) +1
    x_centre <- ((x_range[1]-1) + i) # the left hand edge of the grid cell
    y_centre <- ((y_range[1]-1) + j) # the bottom edge of the grid cell
    grid_centres$centre[count] <- paste0("x", x_centre, "y", y_centre)
    grid_centres$x_grid[count] <- x_centre
    grid_centres$y_grid[count] <- y_centre
    # now calculate lat long
    grid_centres$latitude[count] <- lat[x_centre, y_centre]
    grid_centres$longitude[count] <- lon[x_centre, y_centre]
  }
}

# export as CSV:
write.csv(grid_centres,"/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/grid_centres.csv")

#####

### Create list of cells to be included for anaylsis:
#####
### Cell IDs were determined by masking out the land in QGIS based on the above created grid

# Manual input of cells
qhi_cells <- data.frame(cell_id = rep(NA, 153), x_grid = rep(0, 153), y_grid = rep(0, 153), lat = rep(0, 153), long = rep(0, 153), size = rep(0, 153))
qhi_cells$x_grid <- c(65, rep(66, 2), rep(67, 5), rep(68, 13), rep(69, 14), rep(70, 15), rep(71, 19), rep(72, 21), rep(73, 21), rep(74, 21), rep(75, 21))
qhi_cells$y_grid <- c(231, seq(230, 231, 1), seq(230, 234, 1), seq(226, 238, 1), seq(225, 238, 1), seq(224, 238, 1), seq(220, 238, 1), rep(seq(218, 238, 1), 4))

# Create data frame with cell_id, lat, lon and cell siye collumns
for(i in 1:length(qhi_cells$x_grid)){
  qhi_cells$cell_id[i] <- paste0("x", qhi_cells$x_grid[i], "y", qhi_cells$y_grid[i] )
  qhi_cells$lat[i] <- lat[qhi_cells$x_grid[i], qhi_cells$y_grid[i]]
  qhi_cells$long[i] <- lon[qhi_cells$x_grid[i], qhi_cells$y_grid[i]]
  qhi_cells$size[i] <- pixel_area[qhi_cells$x_grid[i], qhi_cells$y_grid[i]]
}

# Export df as CSV and RDA
write.csv(qhi_cells, "/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_cells.csv")
save(qhi_cells, file = "/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_cells.Rda")
#####

### Now clean up before extracting data across all years
#####
rm(list = ls())
#####

### Extract multi-annual sea-ice extend for the masked area
### !!! WARNING THIS IS LABOUR INTENSE AND CAN TAKE A COUPLE OF HOURS !!!
#####

# Load cell data (mask)
load("/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_cells.Rda")

# Load list of NetCDF files (CDR record downloaded from NSCID)
folder_path <- "/Volumes/csce/biology/users/s1043792/sea_ice/sidads.colorado.edu/pub/DATASETS/NOAA/G02202_v2/north/daily/"
files <- list.files(path = folder_path, recursive = T) # WARNING THIS IS SLOW!
# complete full path to file names in list
files_full <- paste0(folder_path, files)

## Function to calculate sea-ice extend for specified cells from a CDR Sea-Ice NetCDF
# arguments are: file_path and cells_to_extract
# former is an NetCDF file from the CDR and the latter is a data frame with the following collumns at minmum: cell_id, x_grid, y_grid and size 
# akin to the qhi_cells data frame created above

extract_sea_ice_extend <- function(file_path, cells_to_extract) {
  # load file (using ncdf4)
  ncin <- nc_open(file_path)

  # obtain date and create date variables
  ncatt_get(ncin, 0, "time_coverage_start")
  
  file_date <- as.Date(substr(ncatt_get(ncin, 0, "time_coverage_start")$value, 1, 10))
  file_year <- format(file_date, "%Y")
  file_month <- format(file_date, "%b")
  file_day <- format(file_date, "%d")
  file_doy <- format(file_date, "%j")
  
  # load array of sea ice concentrations
  sea_ice_data <- ncvar_get(ncin, "seaice_conc_cdr")
  
  # create empty data frame for sea ice concentration to fill in
  sea_ice_cells <- data.frame(cell_id = cells_to_extract$cell_id, sea_ice_conc = rep(NA, length(cells_to_extract$cell_id)), sea_ice_pa = rep(NA, length(cells_to_extract$cell_id)), sea_ice_extend = rep(NA, length(cells_to_extract$cell_id)))
  # loop through cell data frame and extract sea ice concentrations
  for(i in 1:length(cells_to_extract$cell_id)){
    sea_ice_cells$sea_ice_conc[i] <- sea_ice_data[cells_to_extract$x_grid[i], cells_to_extract$y_grid[i]]
  }

  # determine sea ice presence absence in cell
  sea_ice_cells$sea_ice_pa <- as.numeric(sea_ice_cells$sea_ice_conc >= 0.15)
  # calculate sea ice extend
  sea_ice_cells$sea_ice_extend <- sea_ice_cells$sea_ice_pa * cells_to_extract$size
  # sum up total area
  sea_ice_extend_sum <- sum(sea_ice_cells$sea_ice_extend)
  
  # return all to data frame
  return_df <- data.frame(date = file_date, year = file_year, month = file_month, day = file_day, doy = file_doy, sea_ice_extend = sea_ice_extend_sum)
  
  # and return to common sea_ice_extend data frame
  if (exists("sea_ice_extend")){
    sea_ice_extend <<- rbind (sea_ice_extend, return_df) # NB <<- makes it available globally
  } else {
    sea_ice_extend <<- return_df
  }
  
  # close netcdf file
  nc_close(ncin)

  # return something for good practice
  return(paste0(file_date, " Done."))
}

# run a few tests of the funciton on individual files
extract_sea_ice_extend(files_full[123], qhi_cells)
extract_sea_ice_extend(files_full[250], qhi_cells)
extract_sea_ice_extend(files_full[2534], qhi_cells)
# good, this works well!

# tidy up:
rm(sea_ice_extend)

# now apply to all files in the list (this will take some time)
sapply(files_full, extract_sea_ice_extend, cells_to_extract = qhi_cells)

# write to csv and rdata
write.csv(sea_ice_extend, "/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_sea_ice.csv")
save(sea_ice_extend, file = "/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_sea_ice.Rda")

# For some reason all but the date columns in the sea_ice_extend data frame are saved as factors!
# This is structure is saved in the RData file. :/
load("/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_sea_ice.Rda")
str(sea_ice_extend)
# Let's fix that!

# A quick hack to quickly do the conversion is by loading the CSV isntead
sea_ice_extend <- read.csv("/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_sea_ice.csv")
str(sea_ice_extend)

# get rid of the first column
sea_ice_extend <- sea_ice_extend %>% select(-X)
sea_ice_extend$date <- as.Date(sea_ice_extend$date)

# save the data without the collumns being factors
save(sea_ice_extend, file = "/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_sea_ice.Rda")

### DONE! :)
#####

################################################################################
### Summarise and prepare data for phenology models

### Visualise overall data
#####
load("/Volumes/Users/jakob/Documents/Ecology/PhD/ShrubHub/scripts/users/jassmann/sea_ice/qhi_sea_ice.Rda")

# First check how many data entries are missing
sum(is.na(sea_ice_extend$sea_ice_extend))
length(sea_ice_extend$sea_ice_extend)
# 16 out of a total of 5478, that's not too bad. 
# Data is missing for the following dates:
sea_ice_extend$date[is.na(sea_ice_extend$sea_ice_extend)]

require(ggplot2)
require(dplyr)
require(scales)
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_extend.png", width = 6, height = 4)
#####

### Summarise into annual data by monthly means
#####

# Spring averages (April, May, June) and (May, July)

spring_avg <- sea_ice_extend %>% group_by(year) %>% filter(month == "Apr" | month == "May" | month == "Jun") %>% summarise(amj_mean = mean(sea_ice_extend, na.rm = T)) 
spring_avg_mj <- sea_ice_extend %>% group_by(year) %>% filter(month == "May" | month == "Jun") %>% summarise(mj_mean = mean(sea_ice_extend, na.rm = T)) 
spring_avg$mj_mean <- spring_avg_mj$mj_mean
rm(spring_avg_mj)

## Averages for phenologically important time periods

# salarc -> SALARC earliest snowmelt: doy 130 (10 May), 75% of phenology events: 172 (20 June) with two weeks before: doy 116-171 (26th April - 20 June)
spring_avg_salarc_phen <- sea_ice_extend %>% group_by(year) %>% filter(doy >= 116 & doy <= 171) %>% summarise(salarc_mean = mean(sea_ice_extend, na.rm = T)) 
spring_avg$salarc_mean <- spring_avg_salarc_phen$salarc_mean
rm(spring_avg_salarc_phen)

# erivag -> ERIVAG earlies snowmelt: doy 120 (31 April) 120, 75% of phenology events: doy 152 (1 June) incl. two weeks before: doy 106-152 (16th April - 1 June)

spring_avg_erivag_phen <- sea_ice_extend %>% group_by(year) %>% filter(doy >= 106 & doy <= 152) %>% summarise(erivag_mean = mean(sea_ice_extend, na.rm = T)) 
spring_avg$erivag_mean <- spring_avg_erivag_phen$erivag_mean
rm(spring_avg_erivag_phen)
#####

# Plot average spring extend
#####
library(tidyr)
spring_extend_long <- gather(spring_avg[,1:5], "year", "extend", 2:5)
names(spring_extend_long) <- c("year", "period", "extend")
ggplot(spring_extend_long, aes(x = year , y = extend, group = period, color = period)) + geom_point() + geom_line() +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("\nyear") +
  scale_x_continuous(breaks = seq(2001,2015)) +
  ggtitle("Average spring sea-ice extend by time-period\n") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
rm(spring_extend_long)
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_avg_spring_extend.png", width = 6, height = 4)
#####


################################################################################
### Identify DOY with a given percentage of melt reached

### Identify sea_ice extend at 50% melt, 25% melt, 75%, 80%, 85%, 90% melt
#####
sea_ice_extend_50melt <- sea_ice_extend %>% group_by(year) %>% summarise(min_extend = min(sea_ice_extend, na.rm = T), max_extend = max(sea_ice_extend, na.rm = T), max_doy = max(doy)) %>%
  mutate(sea_50melt_extend = ((max_extend - min_extend)/2 + min_extend), 
         sea_25melt_extend = ((max_extend - min_extend)*0.25 + min_extend), 
         sea_75melt_extend = ((max_extend - min_extend)*0.75 + min_extend),
         sea_80melt_extend = ((max_extend - min_extend)*0.80 + min_extend),
         sea_85melt_extend = ((max_extend - min_extend)*0.85 + min_extend),
         sea_90melt_extend = ((max_extend - min_extend)*0.90 + min_extend),
         sea_95melt_extend = ((max_extend - min_extend)*0.95 + min_extend),
         sea_975melt_extend = ((max_extend - min_extend)*0.975 + min_extend),
         sea_99melt_extend = ((max_extend - min_extend)*0.99 + min_extend)
         )
#####

### Drop below 50% of annual extend (i.e. DOY of drop below 50% melt for that given year)
#####
# identify doy minimum extend is reached, this will be needed as a cut-off later (otherwise we might identify autumn mid-freeze)
# temporary collumn in sea_ice_extend with min anuual sea ice date and then of doy wiht min sea ice
sea_ice_extend$min_ice_extend <- rep(sea_ice_extend_50melt$min_extend, sea_ice_extend_50melt$max_doy)
# temporary df with minimum doy of extend per year
min_extend_doy <- sea_ice_extend %>% filter(sea_ice_extend <= min_ice_extend) %>% group_by(year) %>% summarise(min_extend_doy = min(doy))
# expand into collumn of sea_ice_extend
sea_ice_extend$min_ice_extend_doy <- rep(min_extend_doy$min_extend_doy, sea_ice_extend_50melt$max_doy)
# remove temoporary object
rm(min_extend_doy)

# temporary collumn in sea_ice_extend with 50% sea-ice extend, here we simply expand the sea ice 50% extend value for each year
sea_ice_extend$sea_50melt_extend <- rep(sea_ice_extend_50melt$sea_50melt_extend, sea_ice_extend_50melt$max_doy)

# quality control of the above
sea_ice_extend %>% filter(doy == 1 | doy == 365 | doy == 366) %>% select(year, doy, sea_50melt_extend)
sea_ice_extend %>% group_by(year) %>% summarise(mean_50extend = mean( sea_50melt_extend)) 
# aweseome that worked!

# identify the doy sea ice extend first drops below 50% extend in given year
doy_first_50extend <- sea_ice_extend %>% filter(sea_ice_extend <= sea_50melt_extend) %>% group_by(year) %>% summarise(doy_first_50extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_50extend$date_first_50extend <- rep(as.Date(NA), length(doy_first_50extend$doy_first_50extend))
doy_first_50extend$first_50extend <- rep(as.numeric(NA), length(doy_first_50extend$doy_first_50extend))
for(i in 1:length(doy_first_50extend$doy_first_50extend)){
  first_50extend <- sea_ice_extend %>% filter(year == doy_first_50extend$year[i], doy == doy_first_50extend$doy_first_50extend[i]) %>% select(sea_ice_extend)
  doy_first_50extend$first_50extend[i] <- first_50extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_50extend)
  
  date_first_50extend <- sea_ice_extend %>% filter(year == doy_first_50extend$year[i], doy == doy_first_50extend$doy_first_50extend[i]) %>% select(date)
  doy_first_50extend$date_first_50extend[i] <- date_first_50extend[,1]
  rm(date_first_50extend)
  }
rm(i)

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_50extend, aes(x = date_first_50extend, y = first_50extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below annual 50% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_50melt_ann.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_50extend <- doy_first_50extend$doy_first_50extend 

#####

### Identify dates first drop below a certain absolute extend
### Calculate first doy falling below inter-annual average 50% minimum, 25%, 75%, 80%, 85%, 90%, 95%, 97.5% and 99%
#####
average_50melt <- round(mean(sea_ice_extend_50melt$sea_50melt_extend), 0)
average_25melt <- round(mean(sea_ice_extend_50melt$sea_25melt_extend), 0)
average_75melt <- round(mean(sea_ice_extend_50melt$sea_75melt_extend), 0)
average_80melt <- round(mean(sea_ice_extend_50melt$sea_80melt_extend), 0)
average_85melt <- round(mean(sea_ice_extend_50melt$sea_85melt_extend), 0)
average_90melt <- round(mean(sea_ice_extend_50melt$sea_90melt_extend), 0)
average_95melt <- round(mean(sea_ice_extend_50melt$sea_95melt_extend), 0)
average_975melt <- round(mean(sea_ice_extend_50melt$sea_975melt_extend), 0)
average_99melt <- round(mean(sea_ice_extend_50melt$sea_99melt_extend), 0)
#####

# identify the doy sea ice extend first drops below average 25% extend
######

doy_first_avg_25extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_25melt) %>% group_by(year) %>% summarise(doy_first_avg_25extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_25extend$date_first_avg_25extend <- rep(as.Date(NA), length(doy_first_avg_25extend$doy_first_avg_25extend))
doy_first_avg_25extend$first_avg_25extend <- rep(as.numeric(NA), length(doy_first_avg_25extend$doy_first_avg_25extend))
for(i in 1:length(doy_first_avg_25extend$doy_first_avg_25extend)){
  first_avg_25extend <- sea_ice_extend %>% filter(year == doy_first_avg_25extend$year[i], doy == doy_first_avg_25extend$doy_first_avg_25extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_25extend$first_avg_25extend[i] <- first_avg_25extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_25extend)
  
  date_first_avg_25extend <- sea_ice_extend %>% filter(year == doy_first_avg_25extend$year[i], doy == doy_first_avg_25extend$doy_first_avg_25extend[i]) %>% select(date)
  doy_first_avg_25extend$date_first_avg_25extend[i] <- date_first_avg_25extend[,1]
  rm(date_first_avg_25extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_25extend, aes(x = date_first_avg_25extend, y = first_avg_25extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 25% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_25melt_avg.png", width = 6, height = 4)

# in 2003 the sea ice extend does not drop below the threshold!
length(doy_first_avg_25extend$doy_first_avg_25extend)
# 14 years instead of 15 in the record. 
# Really this renders the 25% threshold unusable!
######

# identify the doy sea ice extend first drops below average 50% extend
######
doy_first_avg_50extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_50melt) %>% group_by(year) %>% summarise(doy_first_avg_50extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_50extend$date_first_avg_50extend <- rep(as.Date(NA), length(doy_first_avg_50extend$doy_first_avg_50extend))
doy_first_avg_50extend$first_avg_50extend <- rep(as.numeric(NA), length(doy_first_avg_50extend$doy_first_avg_50extend))
for(i in 1:length(doy_first_avg_50extend$doy_first_avg_50extend)){
  first_avg_50extend <- sea_ice_extend %>% filter(year == doy_first_avg_50extend$year[i], doy == doy_first_avg_50extend$doy_first_avg_50extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_50extend$first_avg_50extend[i] <- first_avg_50extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_50extend)
  
  date_first_avg_50extend <- sea_ice_extend %>% filter(year == doy_first_avg_50extend$year[i], doy == doy_first_avg_50extend$doy_first_avg_50extend[i]) %>% select(date)
  doy_first_avg_50extend$date_first_avg_50extend[i] <- date_first_avg_50extend[,1]
  rm(date_first_avg_50extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_50extend, aes(x = date_first_avg_50extend, y = first_avg_50extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 50% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_50melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_50extend <- doy_first_avg_50extend$doy_first_avg_50extend 
######

# identify the doy sea ice extend first drops below average 75% extend
######
doy_first_avg_75extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_75melt) %>% group_by(year) %>% summarise(doy_first_avg_75extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_75extend$date_first_avg_75extend <- rep(as.Date(NA), length(doy_first_avg_75extend$doy_first_avg_75extend))
doy_first_avg_75extend$first_avg_75extend <- rep(as.numeric(NA), length(doy_first_avg_75extend$doy_first_avg_75extend))
for(i in 1:length(doy_first_avg_75extend$doy_first_avg_75extend)){
  first_avg_75extend <- sea_ice_extend %>% filter(year == doy_first_avg_75extend$year[i], doy == doy_first_avg_75extend$doy_first_avg_75extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_75extend$first_avg_75extend[i] <- first_avg_75extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_75extend)
  
  date_first_avg_75extend <- sea_ice_extend %>% filter(year == doy_first_avg_75extend$year[i], doy == doy_first_avg_75extend$doy_first_avg_75extend[i]) %>% select(date)
  doy_first_avg_75extend$date_first_avg_75extend[i] <- date_first_avg_75extend[,1]
  rm(date_first_avg_75extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_75extend, aes(x = date_first_avg_75extend, y = first_avg_75extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 75% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_75melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_75extend <- doy_first_avg_75extend$doy_first_avg_75extend 

######

# identify the doy sea ice extend first drops below average 80% extend
######
doy_first_avg_80extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_80melt) %>% group_by(year) %>% summarise(doy_first_avg_80extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_80extend$date_first_avg_80extend <- rep(as.Date(NA), length(doy_first_avg_80extend$doy_first_avg_80extend))
doy_first_avg_80extend$first_avg_80extend <- rep(as.numeric(NA), length(doy_first_avg_80extend$doy_first_avg_80extend))
for(i in 1:length(doy_first_avg_80extend$doy_first_avg_80extend)){
  first_avg_80extend <- sea_ice_extend %>% filter(year == doy_first_avg_80extend$year[i], doy == doy_first_avg_80extend$doy_first_avg_80extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_80extend$first_avg_80extend[i] <- first_avg_80extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_80extend)
  
  date_first_avg_80extend <- sea_ice_extend %>% filter(year == doy_first_avg_80extend$year[i], doy == doy_first_avg_80extend$doy_first_avg_80extend[i]) %>% select(date)
  doy_first_avg_80extend$date_first_avg_80extend[i] <- date_first_avg_80extend[,1]
  rm(date_first_avg_80extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_80extend, aes(x = date_first_avg_80extend, y = first_avg_80extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 80% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_80melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_80extend <- doy_first_avg_80extend$doy_first_avg_80extend 

######

# identify the doy sea ice extend first drops below average 85% extend
######
doy_first_avg_85extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_85melt) %>% group_by(year) %>% summarise(doy_first_avg_85extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_85extend$date_first_avg_85extend <- rep(as.Date(NA), length(doy_first_avg_85extend$doy_first_avg_85extend))
doy_first_avg_85extend$first_avg_85extend <- rep(as.numeric(NA), length(doy_first_avg_85extend$doy_first_avg_85extend))
for(i in 1:length(doy_first_avg_85extend$doy_first_avg_85extend)){
  first_avg_85extend <- sea_ice_extend %>% filter(year == doy_first_avg_85extend$year[i], doy == doy_first_avg_85extend$doy_first_avg_85extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_85extend$first_avg_85extend[i] <- first_avg_85extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_85extend)
  
  date_first_avg_85extend <- sea_ice_extend %>% filter(year == doy_first_avg_85extend$year[i], doy == doy_first_avg_85extend$doy_first_avg_85extend[i]) %>% select(date)
  doy_first_avg_85extend$date_first_avg_85extend[i] <- date_first_avg_85extend[,1]
  rm(date_first_avg_85extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_85extend, aes(x = date_first_avg_85extend, y = first_avg_85extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 85% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_85melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_85extend <- doy_first_avg_85extend$doy_first_avg_85extend 

######

# identify the doy sea ice extend first drops below average 90% extend
######
doy_first_avg_90extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_90melt) %>% group_by(year) %>% summarise(doy_first_avg_90extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_90extend$date_first_avg_90extend <- rep(as.Date(NA), length(doy_first_avg_90extend$doy_first_avg_90extend))
doy_first_avg_90extend$first_avg_90extend <- rep(as.numeric(NA), length(doy_first_avg_90extend$doy_first_avg_90extend))
for(i in 1:length(doy_first_avg_90extend$doy_first_avg_90extend)){
  first_avg_90extend <- sea_ice_extend %>% filter(year == doy_first_avg_90extend$year[i], doy == doy_first_avg_90extend$doy_first_avg_90extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_90extend$first_avg_90extend[i] <- first_avg_90extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_90extend)
  
  date_first_avg_90extend <- sea_ice_extend %>% filter(year == doy_first_avg_90extend$year[i], doy == doy_first_avg_90extend$doy_first_avg_90extend[i]) %>% select(date)
  doy_first_avg_90extend$date_first_avg_90extend[i] <- date_first_avg_90extend[,1]
  rm(date_first_avg_90extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_90extend, aes(x = date_first_avg_90extend, y = first_avg_90extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 90% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_90melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_90extend <- doy_first_avg_90extend$doy_first_avg_90extend 

######

# identify the doy sea ice extend first drops below average 95% extend
######
doy_first_avg_95extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_95melt) %>% group_by(year) %>% summarise(doy_first_avg_95extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_95extend$date_first_avg_95extend <- rep(as.Date(NA), length(doy_first_avg_95extend$doy_first_avg_95extend))
doy_first_avg_95extend$first_avg_95extend <- rep(as.numeric(NA), length(doy_first_avg_95extend$doy_first_avg_95extend))
for(i in 1:length(doy_first_avg_95extend$doy_first_avg_95extend)){
  first_avg_95extend <- sea_ice_extend %>% filter(year == doy_first_avg_95extend$year[i], doy == doy_first_avg_95extend$doy_first_avg_95extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_95extend$first_avg_95extend[i] <- first_avg_95extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_95extend)
  
  date_first_avg_95extend <- sea_ice_extend %>% filter(year == doy_first_avg_95extend$year[i], doy == doy_first_avg_95extend$doy_first_avg_95extend[i]) %>% select(date)
  doy_first_avg_95extend$date_first_avg_95extend[i] <- date_first_avg_95extend[,1]
  rm(date_first_avg_95extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_95extend, aes(x = date_first_avg_95extend, y = first_avg_95extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 95% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_95melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_95extend <- doy_first_avg_95extend$doy_first_avg_95extend 

######

# identify the doy sea ice extend first drops below average 975% extend
######
doy_first_avg_975extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_975melt) %>% group_by(year) %>% summarise(doy_first_avg_975extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_975extend$date_first_avg_975extend <- rep(as.Date(NA), length(doy_first_avg_975extend$doy_first_avg_975extend))
doy_first_avg_975extend$first_avg_975extend <- rep(as.numeric(NA), length(doy_first_avg_975extend$doy_first_avg_975extend))
for(i in 1:length(doy_first_avg_975extend$doy_first_avg_975extend)){
  first_avg_975extend <- sea_ice_extend %>% filter(year == doy_first_avg_975extend$year[i], doy == doy_first_avg_975extend$doy_first_avg_975extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_975extend$first_avg_975extend[i] <- first_avg_975extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_975extend)
  
  date_first_avg_975extend <- sea_ice_extend %>% filter(year == doy_first_avg_975extend$year[i], doy == doy_first_avg_975extend$doy_first_avg_975extend[i]) %>% select(date)
  doy_first_avg_975extend$date_first_avg_975extend[i] <- date_first_avg_975extend[,1]
  rm(date_first_avg_975extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_975extend, aes(x = date_first_avg_975extend, y = first_avg_975extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 97.5% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_975melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_975extend <- doy_first_avg_975extend$doy_first_avg_975extend 

######

# identify the doy sea ice extend first drops below average 99% extend
######
doy_first_avg_99extend <- sea_ice_extend %>% filter(sea_ice_extend <= average_99melt) %>% group_by(year) %>% summarise(doy_first_avg_99extend = min(doy))
# look up sea_ice_extend  and date at above identified doy
# prep columns in df (this is needed for date conversions)
doy_first_avg_99extend$date_first_avg_99extend <- rep(as.Date(NA), length(doy_first_avg_99extend$doy_first_avg_99extend))
doy_first_avg_99extend$first_avg_99extend <- rep(as.numeric(NA), length(doy_first_avg_99extend$doy_first_avg_99extend))
for(i in 1:length(doy_first_avg_99extend$doy_first_avg_99extend)){
  first_avg_99extend <- sea_ice_extend %>% filter(year == doy_first_avg_99extend$year[i], doy == doy_first_avg_99extend$doy_first_avg_99extend[i]) %>% select(sea_ice_extend)
  doy_first_avg_99extend$first_avg_99extend[i] <- first_avg_99extend[,1] # rather clumsy way is needed to convert df into vector
  rm(first_avg_99extend)
  
  date_first_avg_99extend <- sea_ice_extend %>% filter(year == doy_first_avg_99extend$year[i], doy == doy_first_avg_99extend$doy_first_avg_99extend[i]) %>% select(date)
  doy_first_avg_99extend$date_first_avg_99extend[i] <- date_first_avg_99extend[,1]
  rm(date_first_avg_99extend)
}

# quality control by plotting dates on inter-annual sea-ice curve
ggplot(sea_ice_extend, aes(x = date, y = sea_ice_extend)) + geom_line() +
  geom_point(data = doy_first_avg_99extend, aes(x = date_first_avg_99extend, y = first_avg_99extend)) +
  scale_x_date(labels = date_format("%Y"), breaks = "1 year") +
  ylab("Sea-Ice Extend (km2)\n") +
  xlab("") +
  ggtitle("Drop below multi-annual average 99% extend") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_99melt_avg.png", width = 6, height = 4)

# copy into sprin_avg df
spring_avg$doy_first_avg_99extend <- doy_first_avg_99extend$doy_first_avg_99extend 

######

### Clean up section above
#####
rm(doy_first_50extend, sea_ice_extend_50melt, doy_first_avg_25extend, doy_first_avg_50extend, doy_first_avg_75extend, doy_first_avg_80extend, doy_first_avg_85extend, doy_first_avg_90extend, doy_first_avg_95extend, doy_first_avg_975extend, doy_first_avg_99extend)
rm(average_25melt, average_50melt, average_75melt, average_80melt, average_85melt, average_90melt, average_95melt, average_975melt, average_99melt, i)
sea_ice_extend$min_ice_extend <- NULL
sea_ice_extend$min_ice_extend_doy <- NULL
sea_ice_extend$sea_50melt_extend <- NULL
#####

# Plot DOY of melt
#####
library(tidyr)
spring_melt_long <- gather(spring_avg[,c(T, rep(F,4), rep(T,9))], "year", "extend", 2:10)
names(spring_melt_long) <- c("year", "melt", "doy")
ggplot(spring_melt_long, aes(x = year , y = doy, group = melt, color = melt)) + geom_point() + geom_line() +
  ylab("doy\n") +
  xlab("\nyear") +
  scale_x_continuous(breaks = seq(2001,2015)) +
  ggtitle("Day-of-year sea-ice ice melt by extend of melt\n") +
  geom_hline(yintercept = 172, colour = "red") +
  geom_hline(yintercept = 152, colour = "blue") +
  annotate("text", x = 2013, y = 75, colour = "red", label = "75% SALARC") +
  annotate("text", x = 2013, y = 65, colour = "blue", label = "75% ERIVAG") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/qhi_sea_ice_melt_doy.png", width = 6, height = 4)
rm(spring_melt_long)

#####

################################################################################
### Data screening, quality control and export
#####

# Quick correlation matrix of all spring_avg variables
dim(spring_avg)
cormatrix <- cor(spring_avg[,2:dim(spring_avg)[2]])
cormatrix <- rbind(cormatrix, apply(cormatrix, 2, mean)) # calculate mean 
row.names(cormatrix)[dim(cormatrix)[1]] <- "Mean" # give row with means a sensible name 
write.csv(cormatrix, "scripts/users/jassmann/sea_ice/sea_ice_var_cormat.csv") # export
# Drop below 85% correlates the highest with all the other variables, average correlation of 0.85 => Onset of spring melt?!

# Let's look at the inter-annual trend, is there an advance in the 85% onset of melt date? -> Linear model
mod <- lm(spring_avg$year ~ spring_avg$doy_first_avg_85extend)
summary(mod)
# Yes, seems like there is a significant advance in onset of regional melt.

# tidy up
rm(mod)
#####

write.csv(spring_avg, "scripts/users/jassmann/sea_ice/sea_ice_avgs.csv")
save(spring_avg, file = "scripts/users/jassmann/sea_ice/sea_ice_avgs.Rda")

################################################################################
### Corellation with humidity data from Barter Island

### Perpare humidity data
#####

# load barter island data
load("scripts/users/jassmann/sea_ice/barter_island_weather/barter_weather_clean.Rda")

require(dplyr)

# calculate average humidities
# Spring averages (April, May, June) and (May, July)
spring_hum <- barter_weather %>% group_by(year) %>% filter(month == "Apr" | month == "May" | month == "Jun") %>% summarise(amj_hum = mean(hum_avg, na.rm = T)) 
spring_hum_mj <- barter_weather %>% group_by(year) %>% filter(month == "May" | month == "Jun") %>% summarise(mj_hum = mean(hum_avg, na.rm = T)) 
spring_hum$mj_hum <- spring_hum_mj$mj_hum
rm(spring_hum_mj)

# check for NAs in amj
barter_weather %>% group_by(year) %>% filter(month == "Apr" | month == "May" | month == "Jun") %>% summarise(amj_NAs = sum(is.na(hum_avg)), total = sum(is.na(hum_avg) + !is.na(hum_avg))) 
# 2013 and 2015 are missing more than 30% chuck them out
spring_hum$amj_hum[c(13,15)] <- NA 

# check for NAs in mj
barter_weather %>% group_by(year) %>% filter(month == "May" | month == "Jun") %>% summarise(amj_NAs = sum(is.na(hum_avg)), total = sum(is.na(hum_avg) + !is.na(hum_avg))) 
# 2013 and 2015 are missing more 50%+ chuck them out
spring_hum$mj_hum[c(13,15)] <- NA 


## Averages for phenologically important time periods

# salarc -> SALARC earliest snowmelt: doy 130 (10 May), 75% of phenology events: 172 (20 June) with two weeks before: doy 116-171 (26th April - 20 June)
spring_hum_salarc_phen <- barter_weather %>% group_by(year) %>% filter(doy >= 116 & doy <= 171) %>% summarise(salarc_hum = mean(hum_avg, na.rm = T)) 
spring_hum$salarc_hum <- spring_hum_salarc_phen$salarc_hum
rm(spring_hum_salarc_phen)

# erivag -> ERIVAG earlies snowmelt: doy 120 (31 April) 120, 75% of phenology events: doy 152 (1 June) incl. two weeks before: doy 106-152 (16th April - 1 June)

spring_hum_erivag_phen <- barter_weather %>% group_by(year) %>% filter(doy >= 106 & doy <= 152) %>% summarise(erivag_hum = mean(hum_avg, na.rm = T)) 
spring_hum$erivag_hum <- spring_hum_erivag_phen$erivag_hum
rm(spring_hum_erivag_phen)

# check for NAs in salarc_hum
barter_weather %>% group_by(year) %>% filter(doy >= 116 & doy <= 171) %>% summarise(amj_NAs = sum(is.na(hum_avg)), total = sum(is.na(hum_avg) + !is.na(hum_avg))) 
# 10% missing in 2014 not too bad

# check for NAs in erivag_hum
barter_weather %>% group_by(year) %>% filter(doy >= 106 & doy <= 152) %>% summarise(amj_NAs = sum(is.na(hum_avg)), total = sum(is.na(hum_avg) + !is.na(hum_avg))) 
# all good.

# add sea ice melt to humidity df
spring_hum$doy_first_avg_85extend <- c(spring_avg$doy_first_avg_85extend, NA)
spring_hum$doy_first_avg_95extend <- c(spring_avg$doy_first_avg_95extend, NA)
spring_hum$amj_ice_extend <- c(spring_avg$amj_mean, NA)
spring_hum$mj_ice_extend <- c(spring_avg$mj_mean, NA)

# turn year column into numeric
spring_hum$year <- as.numeric(spring_hum$year)
#####

# Plot Average Humidity
#####
require(tidyr)
spring_hum_long <- gather(spring_hum[, 1:5], "year", "avg_hum", 2:5)
names(spring_hum_long) <- c("year", "period", "avg_hum")
ggplot(spring_hum_long, aes(x = year , y = avg_hum, group = period, color = period)) + geom_point() + geom_line() +
  ylab("rel. humidity\n") +
  xlab("\nyear") +
  scale_x_continuous(breaks = seq(2001,2015)) +
  ggtitle("Barter Island average humidity by period in spring\n") +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("scripts/users/jassmann/sea_ice/bart_hum.png", width = 6, height = 4)
rm(spring_hum_long)

### Test correlations
#####
hum_cor_matrix <- cor(spring_hum[,2:dim(spring_hum)[2]], use = "complete.obs")
# 50% correlation between sea ice melt and humidity, not as high as expected

### Write hum_cor_matrix into file
write.csv(hum_cor_matrix, "scripts/users/jassmann/sea_ice/hum_cor_mat.csv")
#####

# EOF