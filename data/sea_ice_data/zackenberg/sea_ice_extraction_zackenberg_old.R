### The purpose of this script is to extract Short Microwave Radiaiton Sea-Ice data from
# the NSCID Cilmate Recrod Data Set and prepare the data for inclusion in the phenology models
# Jakob Assmann (j.assmann@ed.ac.uk) August 2017

# load (some of the) dependencies
library(chron)
library(lattice)
library(ncdf4)
library(dplyr)
library(ggplot2)

# Housekeeping: set folder paths and meta_variables
script_path <- "scripts/users/jassmann/phenology/sea_ice_data/zackenberg/"
site_name <- "zackenberg"

# Site Location
site_lat_lon <- c(74.47, -20.56) # Zackenberg ITEX site

# Polar Stereographic Grid Pixel Area File (binary)
pixel_area_file <- "scripts/users/jassmann/phenology/sea_ice_data/psn25area_v3.dat"

# NSIDC CDR folder path
# folder_path <- "/Volumes/csce/biology/users/s1043792/sea_ice/nsidc_cdr_v3/"
folder_path <- "M:/sea_ice/nsidc_cdr_v3/"

# NB Site Mask needs to be specified and edited manually (Line 200)

################################################################################
### Extract Sea Ice Data from NSCID CDR

### Load sample Climate Data Record NetCDF file, check dimensions, properties and familiarise
#####

# Sample file file path
file_path <- paste0(folder_path, "2017/seaice_conc_daily_nh_f17_20170110_v03r01.nc")

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
# The cell containing Zackenberg is: x 183 and y 297 !
###

# Confirm this by checking lat long ouptut
centre_lat_lon <- c(lat[centre_index], lon[centre_index])
centre_lat_lon

# Clauclate boundary box in xy grid that has the barrow cell in centre (going 10 rows/collumns in each direction)
x_range <- c(centre_index[1] - 10, centre_index[1] + 10)
y_range <- c(centre_index[2] - 10, centre_index[2] + 10)

#####

### Load pixel area data from binary file. 
#####
# prep matrix for storing data
pixel_area <- matrix(nrow = 304, ncol = 448)
# read in binary pixel area file (file path specified in header of this script)
to_read <- file(pixel_area_file, "rb") 
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
write.csv(grid_half_x, paste0(script_path, "/grid_half_", site_name, "_x.csv"))
write.csv(grid_half_y, paste0(script_path, "/grid_half_", site_name, "_y.csv"))

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
write.csv(grid_centres,paste0(script_path, "/grid_centres_", site_name, ".csv"))

# Manual creation of mask required!
# -> improt grid croners (grid halfs x and y) into QHIS 
# then join using Points2One Plugin to create lines features, export as KML files
# then import into google maps to create list of cells to be included by identifiyng 
# all cells that are that are not land or one removed from land. 
# Can also be done in QGIS if a coastal shapefile is available.

#####

### Create list of cells to be included for anaylsis:
#####
### Cell IDs were determined by masking out the land in QGIS / Google Earth 
### as descirbed above

# Manual input of cells that are not masked, x then y low to high 
zackenberg_cells <- data.frame(cell_id = rep(NA, 21 * 6 + 20 + (8+ 4) + 4 + 3 + 1), x_grid = rep(0, 166), y_grid = rep(0, 166), lat = rep(0, 166), long = rep(0, 166), size = rep(0, 166))
zackenberg_cells$x_grid <- c(rep(183, 1), rep(184, 3), rep(185, 4), rep(186, 12), rep(187, 20), rep(188, 21), rep(189, 21), rep(190, 21), rep(191, 21), rep(192, 21), rep(193, 21))
zackenberg_cells$y_grid <- c(287, seq(287, 289, 1), seq(287, 290, 1), seq(287, 290, 1), seq(296, 303, 1), seq(287, 306, 1), rep(seq(287, 307, 1), 6))

# Create fucntion to return data frame with cell_id, lat, lon and cell size collumns
fill_cells_df <- function(cells_df){
  for(i in 1:length(cells_df$x_grid)){
    cells_df$cell_id[i] <- paste0("x", cells_df$x_grid[i], "y", cells_df$y_grid[i] )
    cells_df$lat[i] <- lat[cells_df$x_grid[i], cells_df$y_grid[i]]
    cells_df$long[i] <- lon[cells_df$x_grid[i], cells_df$y_grid[i]]
    cells_df$size[i] <- pixel_area[cells_df$x_grid[i], cells_df$y_grid[i]]
  }
  return(cells_df)
}

# execute function to fill dataframe
zackenberg_cells <- fill_cells_df(zackenberg_cells)

# Export df as CSV and RDA
write.csv(zackenberg_cells, paste0(script_path, "zackenberg_cells.csv"))
save(zackenberg_cells, file = paste0(script_path, "zackenberg_cells.Rda"))
#####

### Now clean up before extracting data across all years
#####
rm(list = ls())
#####

### Extract multi-annual sea-ice extent for the masked area
### !!! WARNING THIS IS LABOUR INTENSE AND CAN TAKE A COUPLE OF HOURS !!!
#####

# Load cell data (mask)
load(paste0(script_path, site_name, "_cells.Rda"))

# Load list of NetCDF files (CDR record downloaded from NSCID)
# folder_path specified in header
files <- list.files(path = folder_path, recursive = T) # WARNING THIS IS SLOW!
# complete full path to file names in list
files_full <- paste0(folder_path, files)

## Function to calculate sea-ice extent for specified cells from a CDR Sea-Ice NetCDF
# arguments are: file_path and cells_to_extract
# former is an NetCDF file from the CDR and the latter is a data frame with the following collumns at minmum: cell_id, x_grid, y_grid and size 
# akin to the barrow_cells data frame created above

extract_sea_ice_extent <- function(file_path, cells_to_extract) {
  # load file (using ncdf4)
  ncin <- nc_open(file_path)

  # obtain date and create date variables
  ncatt_get(ncin, 0, "time_coverage_start")
  
  file_date <- as.Date(substr(ncatt_get(ncin, 0, "time_coverage_start")$value, 1, 10))
  file_year <- format(file_date, "%Y")
  file_month <- format(file_date, "%m")
  file_day <- format(file_date, "%d")
  file_doy <- format(file_date, "%j")
  
  # load array of sea ice concentrations
  sea_ice_data <- ncvar_get(ncin, "seaice_conc_cdr")
  
  # create empty data frame for sea ice concentration to fill in
  sea_ice_cells <- data.frame(cell_id = cells_to_extract$cell_id, sea_ice_conc = rep(NA, length(cells_to_extract$cell_id)), sea_ice_pa = rep(NA, length(cells_to_extract$cell_id)), sea_ice_extent = rep(NA, length(cells_to_extract$cell_id)))
  # loop through cell data frame and extract sea ice concentrations
  for(i in 1:length(cells_to_extract$cell_id)){
    sea_ice_cells$sea_ice_conc[i] <- sea_ice_data[cells_to_extract$x_grid[i], cells_to_extract$y_grid[i]]
  }

  # determine sea ice presence absence in cell
  sea_ice_cells$sea_ice_pa <- as.numeric(sea_ice_cells$sea_ice_conc >= 0.15)
  # calculate sea ice extent
  sea_ice_cells$sea_ice_extent <- sea_ice_cells$sea_ice_pa * cells_to_extract$size
  # sum up total area and round to full number
  sea_ice_extent_sum <- round(sum(sea_ice_cells$sea_ice_extent))
  
  # return all to data frame
  return_df <- data.frame(date = file_date, year = file_year, month = file_month, day = file_day, doy = file_doy, sea_ice_extent = sea_ice_extent_sum)
  
  # and return to common sea_ice_extent data frame
  if (exists("sea_ice_extent")){
    sea_ice_extent <<- rbind (sea_ice_extent, return_df) # NB <<- makes it available globally
  } else {
    sea_ice_extent <<- return_df
  }
  
  # close netcdf file
  nc_close(ncin)

  # print status update
  cat(as.character(file_date), "")
  
  # return something for good practice
  return(paste0(file_date, " Done."))
}

# run a few tests of the funciton on individual files
extract_sea_ice_extent(files_full[123], get(paste0(site_name, "_cells")))
extract_sea_ice_extent(files_full[250], get(paste0(site_name, "_cells")))
extract_sea_ice_extent(files_full[2534], get(paste0(site_name, "_cells")))
# good, this works well!

# tidy up:
rm(sea_ice_extent)

# now apply to all files in the list (this will take some time)
sapply(files_full, extract_sea_ice_extent, cells_to_extract = get(paste0(site_name, "_cells")))

# add collumn with site_name
sea_ice_extent$site_name <- toupper(site_name)

# re-arrange collumns
sea_ice_extent <- sea_ice_extent[c(7,1:6)]

# write to csv and rdata
write.csv(sea_ice_extent, paste0(script_path, site_name, "_sea_ice_extent.csv"))
save(sea_ice_extent, file = paste0(script_path, site_name, "_sea_ice_extent.Rda"))

# For some reason all but the date columns in the sea_ice_extent data frame are saved as factors!
# This is structure is saved in the RData file. :/
load(paste0(script_path, site_name, "_sea_ice_extent.Rda"))
str(sea_ice_extent)
# Let's fix that!

# A quick hack to quickly do the conversion is by loading the CSV isntead
sea_ice_extent <- read.csv(paste0(script_path, site_name, "_sea_ice_extent.csv"))
str(sea_ice_extent)

# get rid of the first column
sea_ice_extent <- sea_ice_extent %>% select(-X)
sea_ice_extent$date <- as.Date(sea_ice_extent$date)

# save the data without the collumns being factors
save(sea_ice_extent, file = paste0(script_path, site_name, "_sea_ice_extent.Rda"))

### DONE! :)
#####

### EOF