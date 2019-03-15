### Preparation of QHI phenology data for the Coastal Phenology Analysis
## 8 Feb 2018 Jakob Assmann j.assmann@ed.ac.uk

# dependencies
library(dplyr)

# Housekeeping
script_path <- "phenology_data/qhi/"

qhi_phen <- read.csv("phenology_data/qhi/qiki_phen_with_before_2017.csv")

### The aim is to create a df with the following structure and 
### only containing: ERIVAG flowering, SALARC green_up and DRYINT flowering
### 
# * site_name (SITE)
# * spp (SPECIES)
# * phen_stage (flowering, green_up)
# * site_spp_phase (SITE_SPECIES_phen)
# * spp_phase (SPECIES_phen)
# * plot_id (plot_id)
# * year (YYYY)
# * doy (J)
# * prior_obs (J)
# * snowmelt (J)

### ERIVAG flowering
erivag_flowering <- qhi_phen %>% filter(Spp == "ERIVAG") %>% 
  select (Spp, Plot.ID, Year, P2, P2_before, P1)
erivag_flowering$phen_stage <- "flowering"  
names(erivag_flowering) <- c("spp", "plot_id", "year", "doy", "prior_visit", "snowmelt", "phen_stage")
erivag_flowering$site_name <- "QHI"
erivag_flowering <- erivag_flowering %>% 
  mutate(spp_phen = paste0(spp, "_", phen_stage), site_spp_phen = paste0(site_name, "_", spp, "_", phen_stage) )
erivag_flowering <- erivag_flowering[, c(8,1,7,10,9,2,3,4,5,6)]

### SALARC green_up
salarc_green_up <- qhi_phen %>% filter(Spp == "SALARC") %>% 
  select (Spp, Plot.ID, Year, P2, P2_before, P1)
salarc_green_up$phen_stage <- "green_up"  
names(salarc_green_up) <- c("spp", "plot_id", "year", "doy", "prior_visit", "snowmelt", "phen_stage")
salarc_green_up$site_name <- "QHI"
salarc_green_up <- salarc_green_up %>% 
  mutate(spp_phen = paste0(spp, "_", phen_stage), site_spp_phen = paste0(site_name, "_", spp, "_", phen_stage) )
salarc_green_up <- salarc_green_up[, c(8,1,7,10,9,2,3,4,5,6)]

### DRYINT flowering
dryint_flowering <- qhi_phen %>% filter(Spp == "DRYINT") %>% 
  select (Spp, Plot.ID, Year, P2, P2_before, P1)
dryint_flowering$phen_stage <- "flowering"  
names(dryint_flowering) <- c("spp", "plot_id", "year", "doy", "prior_visit", "snowmelt", "phen_stage")
dryint_flowering$site_name <- "QHI"
dryint_flowering <- dryint_flowering %>% 
  mutate(spp_phen = paste0(spp, "_", phen_stage), site_spp_phen = paste0(site_name, "_", spp, "_", phen_stage) )
dryint_flowering <- dryint_flowering[, c(8,1,7,10,9,2,3,4,5,6)]

### Check consistency 
names(erivag_flowering) == names(salarc_green_up)
names(salarc_green_up) == names(dryint_flowering)

### Combine into one dataframe
qhi_phen <- rbind(erivag_flowering, salarc_green_up, dryint_flowering)

### Create unique plot identifiers
qhi_phen$plot_id <- paste0(qhi_phen$spp, "_", qhi_phen$plot_id)

### Check whether there are intevals wider than 10
sum(qhi_phen$doy - qhi_phen$prior_visit > 10, na.rm = T)
qhi_phen[which(qhi_phen$doy - qhi_phen$prior_visit > 10),]
mean(qhi_phen[which(qhi_phen$doy - qhi_phen$prior_visit > 10),]$doy -
       qhi_phen[which(qhi_phen$doy - qhi_phen$prior_visit > 10),]$prior_visit)
# 30 day mean window for these two data points, ach well...
# wonder what happened that year?

### Export data
write.csv(qhi_phen, paste0(script_path, "qhi_phen.csv"), row.names=FALSE)
