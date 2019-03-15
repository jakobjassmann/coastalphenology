### Preparation of zackenberg phenology data for the Coastal Phenology Analysis
## 8 Feb 2018 Jakob Assmann j.assmann@ed.ac.uk

# Blarb for website As part of Isla Myers-Smith's group we are currently working on an attribution analysis of the environmental factors that explain variation in the phenology observations at coastal sites in the ITEX control dataset. My plan is to update the Zackenberg phenology data with the most recent observations.

# dependencies
library(dplyr)

# Housekeeping
script_path <- "phenology_data/zackenberg/"

zackenberg_phen <- read.csv("phenology_data/zackenberg/Zackenberg_plotleveldata_1994-2014.csv")
### The aim is to create a df with the following structure and 
### only containing: CASTET green_up  DUPFIS green_up LUZARC flowering "LUZARC green_up
### POAARC green_up  SALROT flowering SALROT green_up
### 
# * site_name (SITE)
# * spp (SPECIES)
# * phen_stage (flowering, green_up)
# * site_spp_phen (SITE_SPECIES_phen)
# * spp_phen (SPECIES_phen)
# * plot_id (plot_id)
# * year (YYYY)
# * doy (J)
# * prior_obs (J)
# * snowmelt (J)
spp_phase_of_interest <- c("CASTET_flowering",
                           "DRYOCT_flowering",
                           "PAPRAD_flowering",
                           "SALARC_flowering",
                           "SAXOPP_flowering",
                           "SILACA_flowering")

# Filter data set to contain only spp_phase combinations of interest
zackenberg_phen$phenophase <- as.character(zackenberg_phen$phenophase)
zackenberg_phen$phenophase[zackenberg_phen$phenophase == "Green"] <- "green_up"
zackenberg_phen$phenophase[zackenberg_phen$phenophase == "Flower"] <- "flowering"
zackenberg_phen <- zackenberg_phen %>% mutate(spp_phen = paste0(spp, "_", phenophase)) %>%
  filter(spp_phen %in% spp_phase_of_interest) %>%
  mutate(site_spp_phen = paste0(site_name, "_", spp, "_", phenophase))

### The data has no snowmelt data!
### Let's grab if from Janet's data
# load Janet's ITEX data for the snow melt data for Zackengberg
itex_phen <- read.csv("scripts/users/jassmann/phenology/phenology_data/ITEX_data_Janet/CCIN12722_20171116_Arctic_phenology_database_1992-2014.csv") %>%
  filter(site_name == "ZACKENBERG")
itex_phen$phenophase <- as.character(itex_phen$phenophase)
itex_phen$phenophase[itex_phen$phenophase == "Green"] <- "green_up"
itex_phen$phenophase[itex_phen$phenophase == "Flower"] <- "flowering"
itex_phen <- itex_phen %>% 
  mutate(spp_phen = paste0(spp, "_", phenophase)) %>%
  filter(spp_phen %in% spp_phase_of_interest)
unique(zackenberg_phen$year)
unique(itex_phen$year)

## We are three years short of snowmelt data only select what's available
zackenberg_phen <- zackenberg_phen %>% filter(year %in% itex_phen$year)

# empty snomwelt collumn
zackenberg_phen$snowmelt <- NA
for (i in min(unique(zackenberg_phen$year)):max(unique(zackenberg_phen$year))) {
  print(paste(i, unique(as.character(itex_phen$snowfree_date[itex_phen$year == i]))))
  zackenberg_phen$snowmelt[zackenberg_phen$year == i] <- unique(as.character(itex_phen$snowfree_date[itex_phen$year == i]))
}

## create Collum with day prior
# use loop for ease of programming, not speed

# create blank collumn in data frame
zackenberg_phen$prior_visit <- NA

# minimum snomwlet on record
zackenberg_phen$snowmelt <- as.numeric(zackenberg_phen$snowmelt)
min_snowmelt <- zackenberg_phen %>% 
  summarise(mean_snowmelt = round(min(snowmelt, na.rm =  T))) %>%
  as.numeric()

for (i in 1:nrow(zackenberg_phen)){
  # extract counter variables
  current_year <- zackenberg_phen$year[i]
  doy <- zackenberg_phen$day[i]
  # filter observations by year, then subset observations prior recorded phenology date
  # get unique ones and slect the maximum value
  # This assumes if one observation was made that year, all plots were visited
  max_pior <- max(unique(zackenberg_phen[zackenberg_phen$year == current_year,]$day[zackenberg_phen[zackenberg_phen$year == current_year,]$day < doy]))
  # if no prior observation was recorded set to 10 days prior
  if(!is.finite(max_pior)) {max_pior <- min_snowmelt}
  # assign to prio_obs collumn
  zackenberg_phen$prior_visit[i] <- max_pior
  #print(max_pior)
}

### Quality Control 
sum(!(zackenberg_phen$prior_obs < zackenberg_phen$day), na.rm = T) # Is 0 (good!)
sum((zackenberg_phen$day - zackenberg_phen$prior_obs) > 10) 
# Brilliant, looking good

# rename collumns
names(zackenberg_phen)[1] <- "site_name"
names(zackenberg_phen)[2] <- "plot_id"
names(zackenberg_phen)[5] <- "phen_stage"
names(zackenberg_phen)[6] <- "doy"

names(zackenberg_phen)

# throw out all unwanted collumns and reorder
zackenberg_phen <- zackenberg_phen %>% select(site_name, spp, phen_stage, site_spp_phen, spp_phen, plot_id, year, doy, prior_visit, snowmelt)

### Export data
write.csv(zackenberg_phen, paste0(script_path, "zackenberg_phen.csv"), row.names=FALSE)
names(zackenberg_phen)
