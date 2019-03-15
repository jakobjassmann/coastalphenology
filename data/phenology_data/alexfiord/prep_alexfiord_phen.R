### Preparation of ALEXFIORD phenology data for the Coastal Phenology Analysis
## 8 Feb 2018 Jakob Assmann j.assmann@ed.ac.uk

# dependencies
library(dplyr)

# Housekeeping
script_path <- "phenology_data/alexfiord/"

alexfiord_phen <- read.csv("../ITEX_data_Janet/CCIN12722_20171116_Arctic_phenology_database_1992-2014.csv") %>%
  filter(site_name == "ALEXFIORD")

### The aim is to create a df with the following structure and 
### the bleow species
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
spp_phase_of_interest <- c("DRYINT_flowering",
                           "DRYINT_green_up",
                           "LUZSPP_flowering",
                           "OXYDIG_flowering",
                           "OXYDIG_green_up",
                           "PAPRAD_flowering",
                           "PAPRAD_green_up",
                           "SALARC_flowering")

# Filter data set to contain only spp_phase combinations of interest
alexfiord_phen$phenophase <- as.character(alexfiord_phen$phenophase)
alexfiord_phen$phenophase[alexfiord_phen$phenophase == "Green"] <- "green_up"
alexfiord_phen$phenophase[alexfiord_phen$phenophase == "Flower"] <- "flowering"
alexfiord_phen <- alexfiord_phen %>% mutate(spp_phen = paste0(spp, "_", phenophase)) %>%
  filter(spp_phen %in% spp_phase_of_interest) %>%
  mutate(site_spp_phen = paste0(site_name, "_", spp, "_", phenophase))

## create Collum with day prior
# use loop for ease of programming, not speed

# create blank collumn in data frame
alexfiord_phen$prior_visit <- NA

# extract min snomwelt across time period of stuy
# for this we need to clean up the snowmelt dates column
alexfiord_phen[alexfiord_phen$snowfree_date == "na" | 
                 alexfiord_phen$snowfree_date == " ",]$snowfree_date <- NA 
alexfiord_phen$snowfree_date <- as.numeric(
  as.character(alexfiord_phen$snowfree_date))
# min snomwlet
min_snowmelt <- alexfiord_phen %>% 
  summarise(mean_snowmelt = round(min(snowfree_date, na.rm =  T))) %>%
  as.numeric()

# fill in prior_visits column by looping through the day vector
# (yes I know it could have been mapped with apply, but I was lazy... :)
for (i in 1:nrow(alexfiord_phen)){
  # extract counter variables
  current_year <- alexfiord_phen$year[i]
  doy <- alexfiord_phen$day[i]
  # filter observations by year, then subset observations prior recorded phenology date
  # get unique ones and slect the maximum value
  # This assumes if one observation was made that year, all plots were visited
  max_pior <- max(unique(alexfiord_phen[alexfiord_phen$year == current_year,]$day[alexfiord_phen[alexfiord_phen$year == current_year,]$day < doy]))
  # if no prior observation was recorded set to mean snowmelt
  if(!is.finite(max_pior)) {max_pior <- min_snowmelt}
  # assign to prio_obs collumn
  alexfiord_phen$prior_visit[i] <- max_pior
  #print(max_pior)
}

### Quality Control 
sum(!(alexfiord_phen$prior_visit < alexfiord_phen$day), na.rm = T) # Is 0 (good!)
sum((alexfiord_phen$day - alexfiord_phen$prior_visit) > 10)
# 556 values where the interval is wider than 10 days
# insepct them in detail
sub_test <- alexfiord_phen[
  which((alexfiord_phen$day - alexfiord_phen$prior_visit) > 10),] %>%
  select(spp_phen, snowfree_date, prior_visit, day)
sub_test <- sub_test[!is.na(sub_test$snowfree_date),]
mean(sub_test$day - sub_test$prior_visit)
# Okay, there are three cases where snomwlet is NA.
# The rest are all beginning of the season estimates with average interval 
# sizes of 17.16 days. This seems good.
rm(sub_test)

# rename collumns
names(alexfiord_phen)[1] <- "site_name"
names(alexfiord_phen)[10] <- "plot_id"
names(alexfiord_phen)[17] <- "phen_stage"
names(alexfiord_phen)[18] <- "doy"
names(alexfiord_phen)[23] <- "snowmelt"
names(alexfiord_phen)

# throw out all unwanted collumns and reorder
alexfiord_phen <- alexfiord_phen %>% 
  select(site_name, spp, phen_stage, site_spp_phen, spp_phen, plot_id, 
         year, doy, prior_visit, snowmelt) %>% 
  filter(year != 1994) # no snowmelt data avialble for 1994!

### Export data
write.csv(alexfiord_phen, paste0(script_path, "alexfiord_phen.csv"), row.names=FALSE)
names(alexfiord_phen)
