# coastal phenology site map
# Jakob Assmann 7 March 2018

# Dependencies
library(maps)
library(raster)
library(maptools)

# Load Site coordinates
site_coordinates <- read.csv("scripts/users/jassmann/phenology/coastal_phenology/site_coordinates.csv")
sites <- SpatialPoints(coords = site_coordinates[,2:3], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))

# Create world map using maps
wrld <- map(plot=FALSE, interior=FALSE, wrap=TRUE, ylim=c(55, 90), xlim=c(-180, 180))
# Transform map to SpatialLines object
wrld_sp <- map2SpatialLines(wrld)
# Set CRS to lat lon (conversion does not carry it across)
proj4string(wrld_sp) <- CRS("+proj=longlat")

# Project map and sites into Canadian Centred polar projection.
laea_wrld_sp <- spTransform(wrld_sp, CRS("+proj=laea +lat_0=90 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
laea_sites <- spTransform(sites,CRS("+proj=laea +lat_0=90 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") )

# Plot map
png(filename = "scripts/users/jassmann/phenology/coastal_phenology/coastal_map.png", width = 1200, height = 1200, res = 600, pointsize = 1)
plot(laea_wrld_sp, col = 'black', lwd = 0.3, cex= 1)
# Plot Sites onto map
plot(laea_sites, pch = 21, cex = 10, 
     col = "black", 
     bg = c("#324D5CFF", 
            "#46B29DFF", 
            "#C2A33EFF",
            "#E37B40FF"),
     add = T)
# Add labels
text(laea_sites, c("Alexandra Fiord",
                   "Utqiaġvik",
                   "Qikiqtaruk",
                   "Zackenberg"), font = 2, cex = 5.7,
     col = c("#324D5CFF", 
             "#46B29DFF", 
             "#C2A33EFF",
             "#E37B40FF"),
     pos = c(4,2,2,4), offset = 3, halo = T, hw = 0.8)
dev.off()

