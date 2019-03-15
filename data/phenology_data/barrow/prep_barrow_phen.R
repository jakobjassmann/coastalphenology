### Preparation of BARROW phenology data for the Coastal Phenology Analysis
## 8 Feb 2018 Jakob Assmann j.assmann@ed.ac.uk

# dependencies
library(dplyr)

# Housekeeping
script_path <- "/phenology_data/barrow/"

barrow_phen <- read.csv(".../ITEX_data_Janet/CCIN12722_20171116_Arctic_phenology_database_1992-2014.csv") %>%
  filter(site_name == "BARROW")

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
spp_phase_of_interest <- c("CASTET_green_up",
                           "DUPFIS_green_up",
                           "LUZARC_flowering",
                           "LUZARC_green_up",
                           "POAARC_green_up",
                           "SALROT_flowering",
                           "SALROT_green_up")

# Filter data set to contain only spp_phase combinations of interest
barrow_phen$phenophase <- as.character(barrow_phen$phenophase)
barrow_phen$phenophase[barrow_phen$phenophase == "Green"] <- "green_up"
barrow_phen$phenophase[barrow_phen$phenophase == "Flower"] <- "flowering"
barrow_phen <- barrow_phen %>% mutate(spp_phen = paste0(spp, "_", phenophase)) %>%
  filter(spp_phen %in% spp_phase_of_interest) %>%
  mutate(site_spp_phen = paste0(site_name, "_", spp, "_", phenophase))

## create Collum with day prior
# use loop for ease of programming, not speed

# create blank collumn in data frame
barrow_phen$prior_visit <- NA

# extract min snomwelt across time period of stuy
barrow_phen$snowfree_date <- as.numeric(
  as.character(barrow_phen$snowfree_date))
# min snomwlet
min_snowmelt <- barrow_phen %>% 
  summarise(mean_snowmelt = round(min(snowfree_date, na.rm =  T))) %>%
  as.numeric()

for (i in 1:nrow(barrow_phen)){
  # extract counter variables
  current_year <- barrow_phen$year[i]
  doy <- barrow_phen$day[i]
  # filter observations by year, then subset observations prior recorded phenology date
  # get unique ones and slect the maximum value
  # This assumes if one observation was made that year, all plots were visited
  max_pior <- max(unique(barrow_phen[barrow_phen$year == current_year,]$day[barrow_phen[barrow_phen$year == current_year,]$day < doy]))
  # if no prior observation was recorded set to 10 days prior
  if(!is.finite(max_pior)) {max_pior <- min_snowmelt}
  # assign to prio_obs collumn
  barrow_phen$prior_visit[i] <- max_pior
  #print(max_pior)
}

### Quality Control 
sum(!(barrow_phen$prior_visit < barrow_phen$day), na.rm = T) # Is 0 (good!)
sum((barrow_phen$day - barrow_phen$prior_visit) > 10) 
# 259 days with intervals wider than 10 days.
sub_test <- barrow_phen[(barrow_phen$day - barrow_phen$prior_visit) > 10,] %>%
  select(spp_phen, snowfree_date, prior_visit, day)
mean(sub_test$day - sub_test$prior_visit)
# mean interval width: 18.44 days. Sounds good!

# rename collumns
names(barrow_phen)[1] <- "site_name"
names(barrow_phen)[10] <- "plot_id"
names(barrow_phen)[17] <- "phen_stage"
names(barrow_phen)[18] <- "doy"
names(barrow_phen)[23] <- "snowmelt"
names(barrow_phen)

# throw out all unwanted collumns and reorder
barrow_phen <- barrow_phen %>% select(site_name, spp, phen_stage, site_spp_phen, spp_phen, plot_id, year, doy, prior_visit, snowmelt)

### Export data
write.csv(barrow_phen, paste0(script_path, "barrow_phen.csv"), row.names=FALSE)
names(barrow_phen)
