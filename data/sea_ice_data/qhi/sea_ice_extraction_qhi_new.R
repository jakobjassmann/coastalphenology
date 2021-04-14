### The purpose of this script is to extract Short Microwave Radiaiton Sea-Ice data from
# the NSCID Cilmate Recrod Data Set and prepare the data for inclusion in the phenology models
# Jakob Assmann (j.assmann@ed.ac.uk) August 2017
# Updated Feb 2018

# load (some of the) dependencies
library(chron)
library(lattice)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(raster)

# house keeping (unix and windows file paths)
script_path <- "data/sea_ice_data/qhi/"
# NSIDC CDR folder of daily observations
# nsidc_cdr_path <- "Volumes/csce/biology/users/s1043792/sea_ice/nsidc_cdr_v3/"
nsidc_cdr_path <- "D:/_TemporaryUserFiles/JakobAssmann/CandiceSeaIce/NSIDC_CDR/daily/"
# pixel area file
pixel_area_file <- "D:/_TemporaryUserFiles/JakobAssmann/CandiceSeaIce/NSIDC_CDR/psn25area_v3.dat"
# site id for files names
site_id <- "qhi"

################################################################################
### Extract Sea Ice Data from NSIDC CDR

### Load sample Climate Data Record NetCDF file and familiarise
#####
# Set a file path to a NetCDF file
file_path <- paste0(nsidc_cdr_path, "2017/seaice_conc_daily_nh_f17_20170110_v03r01.nc")

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
pixel_area_new <- matrix(nrow = 448, ncol = 304)
# read in binary pixel area file
to_read <- file(pixel_area_file, "rb")
# read in pixel area values integer by integer
for (i in 1:448){
  for (j in 1:304){
    pixel_area_new[i, j] <- readBin(to_read, integer(), endian = "little", size = 4) / 1000 
    # NB the "endian" setting is really important, all other files provided by the NSCID with the pixel area file (masks etc) have a little endian, so we can assume that applies to the pixel area file too (though not specified). Also tested big endian to be sure and results do not make sense!
  }
}
# check whether there is anything left (if at end the value will be NULL)
readBin(to_read, integer(), endian = "little", size = 4)
# good!
# close file:
close(to_read)

# check whether all values got read:
pixel_area_new[448,304] # this is the final cell in the data set

plot(raster(pixel_area_new))
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
write.csv(grid_half_x, paste0(script_path, "grid_half_x.csv"))
write.csv(grid_half_y, paste0(script_path, "grid_half_y.csv"))

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
write.csv(grid_centres,paste0(script_path, "grid_centres.csv"))

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
  qhi_cells$size[i] <- pixel_area_new[cells_df$y_grid[i], cells_df$x_grid[i]]
}

# Export df as CSV and RDA
write.csv(qhi_cells,"data/sea_ice_data/qhi/qhi_cells_new.csv")
save(qhi_cells, file = "data/sea_ice_data/qhi/qhi_cells_new.Rda")
#####

### Now clean up before extracting data across all years
#####
rm(list = ls())
#####

### Extract multi-annual sea-ice extend for the masked area
### !!! WARNING THIS IS LABOUR INTENSE AND CAN TAKE A COUPLE OF HOURS !!!
#####

# Load cell data (mask)
load("data/sea_ice_data/qhi/qhi_cells.Rda")

# Load list of NetCDF files (CDR record downloaded from NSCID)
folder_path <- nsidc_cdr_path
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
  file_month <- format(file_date, "%m")
  file_day <- format(file_date, "%d")
  file_doy <- format(file_date, "%j")
  
  # Status
  cat(paste0("Processing file date: ", file_date))
  
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

  # Status
  cat(" Done.")
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

# add collumn with site_name
sea_ice_extend$site_name <- toupper("qhi")


# write to csv and rdata
write.csv(sea_ice_extend, "data/sea_ice_data/qhi/qhi_sea_ice_new.csv")
save(sea_ice_extend, file = "data/sea_ice_data/qhi/qhi_sea_ice_new.Rda")

# For some reason all but the date columns in the sea_ice_extend data frame are saved as factors!
# This is structure is saved in the RData file. :/

# A quick hack to quickly do the conversion is by loading the CSV isntead
sea_ice_extent <- read.csv("data/sea_ice_data/qhi/qhi_sea_ice_new.csv")
str(sea_ice_extent)

# get rid of the first column and correct spelling mistake
sea_ice_extent <- sea_ice_extent %>% dplyr::select(-X)
sea_ice_extent$date <- as.Date(sea_ice_extent$date)
names(sea_ice_extent)[6] <- "sea_ice_extent"

# save the data without the collumns being factors
write.csv(sea_ice_extent, "data/sea_ice_data/qhi/qhi_sea_ice_new.csv", row.names = F)
save(sea_ice_extent, file = "data/sea_ice_data/qhi/qhi_sea_ice_new.Rda")

### DONE! :)
### EOF