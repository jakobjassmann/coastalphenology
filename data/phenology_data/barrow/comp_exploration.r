### Explore community compostion for Barrow

# Dependencies
library(dplyr)

# Load community compositon data
load("data/ITEX_ALL.RData")
names(itex.all)

barrow <- itex.all %>% filter(SITE == "BARROW") %>% 
  group_by(SUBSITE, GENUS, SPECIES) %>%
  summarise(mean_abundance = mean(Abundance)) %>% 
  filter(mean_abundance >= 10)

# Load Janet's ITEX dataset
itex_phen <- read.csv("scripts/users/jassmann/phenology/phenology_data/ITEX_data_Janet/CCIN12722_20171116_Arctic_phenology_database_1992-2014.csv")
sites_of_interest <- c("ALEXFIORD", "BARROW", "ZACKENBERG", "QHI")
phenology <- itex_phen %>% filter(site_name %in% sites_of_interest) %>% group_by(site_name)
site_spp_snow_melt <- phenology %>% group_by(site_name,spp,phenophase, genus, species) %>% summarise(mean_snow_melt = mean(as.numeric(as.character(snowfree_date)), na.rm = T), mean_phenology = mean(day, na.rm = T), difference = mean_phenology - mean_snow_melt, min_year = min(year, na.rm = T), max_year = max(year, na.rm = T), n_years = max_year-min_year+1, n_individ = length(as.character(unique(plot))))
spring_spp <- site_spp_snow_melt %>% filter(difference <= 30 & phenophase != "FlowerEnd" & n_years > 3)

spring_spp_barrow <- spring_spp %>% filter(site_name == "BARROW")

spring_spp_barrow[which(spring_spp_barrow$species %in% barrow$SPECIES),]
barrow[which(barrow$SPECIES %in% spring_spp_barrow$species),]

## Selectes species phenostage combinations are:
selected <- paste0(spring_spp_barrow[which(spring_spp_barrow$species %in% barrow$SPECIES),]$spp, "_",
                   spring_spp_barrow[which(spring_spp_barrow$species %in% barrow$SPECIES),]$phenophase)
selected
