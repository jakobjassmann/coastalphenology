### Coastal Phenology Manuscript Attribution Analysis
# Jakob Assmann, j.assmann@ed.ac.uk 20 Feb 2018

# Dependencies
require(dplyr)
require(ggplot2)
require(MCMCglmm)
require(usdm)
require(snow)
require(doSNOW)

# Housekeeping
# script_path <- "scripts/users/jassmann/phenology/coastal_phenology/"

# Load data
load("data/coastal_phen.Rda")

# set date string
today <- "2021-02-01"

### Quality Control ----

## Check for NAs and quality in predictors and responses

# Phenology observations, upper and lower bounds
unique(coastal_phen$doy)
unique(coastal_phen$prior_visit)
sum(is.na(coastal_phen$doy))
sum(is.na(coastal_phen$prior_visit))
# remove NAs (211 values)
coastal_phen <- coastal_phen[!is.na(coastal_phen$doy),]
# nice!

# Snowmelt dates
unique(coastal_phen$snowmelt)
sum(is.na(coastal_phen$snowmelt))
# Some NAs in there (153)
coastal_phen <- coastal_phen[!is.na(coastal_phen$snowmelt),]
# Nice!

# Remaining predictors
unique(coastal_phen$phase_temp)
unique(coastal_phen$phase_temp_unmodified)
unique(coastal_phen$phase_temp_mean)
unique(coastal_phen$spring_temp)
unique(coastal_phen$apr_temp)
unique(coastal_phen$may_temp)
unique(coastal_phen$jun_temp)
unique(coastal_phen$onset_ice_melt)
unique(coastal_phen$spring_extent)

sum(is.na(coastal_phen$phase_temp))
sum(is.na(coastal_phen$spring_temp))
sum(is.na(coastal_phen$apr_temp))
sum(is.na(coastal_phen$may_temp))
sum(is.na(coastal_phen$jun_temp))
sum(is.na(coastal_phen$onset_ice_melt))
sum(is.na(coastal_phen$spring_extent))
# Okay. All good.

# Remaining years on record
years_on_record <- coastal_phen %>% group_by(site_spp_phen) %>% summarise(min_year = min(year), max_year = max(year))
years_on_record_site <- coastal_phen %>% group_by(site_name) %>% summarise(min_year = min(year), max_year = max(year))

# Create mean snowmelt variable for all sites but Zackenberg
coastal_phen$snowmelt_average <- NA
mean_snowmelt <- coastal_phen %>% group_by(site_name, year) %>% summarise(mean_snowmelt = round(mean(snowmelt),2))
for(year_counter in unique(coastal_phen[coastal_phen$site_name == "ALEXFIORD",]$year)){
  coastal_phen[coastal_phen$site_name == "ALEXFIORD" & coastal_phen$year == year_counter ,]$snowmelt_average <- mean_snowmelt[mean_snowmelt$site_name == "ALEXFIORD" & mean_snowmelt$year == year_counter,]$mean_snowmelt
}
for(year_counter in unique(coastal_phen[coastal_phen$site_name == "BARROW",]$year)){
coastal_phen[coastal_phen$site_name == "BARROW" & coastal_phen$year == year_counter ,]$snowmelt_average <- mean_snowmelt[mean_snowmelt$site_name == "BARROW" & mean_snowmelt$year == year_counter,]$mean_snowmelt
}
for(year_counter in sort(unique(coastal_phen[coastal_phen$site_name == "QHI",]$year))){
  coastal_phen[coastal_phen$site_name == "QHI" & coastal_phen$year == year_counter ,]$snowmelt_average <- mean_snowmelt[mean_snowmelt$site_name == "QHI" & mean_snowmelt$year == year_counter,]$mean_snowmelt
}
for(year_counter in unique(coastal_phen[coastal_phen$site_name == "ZACKENBERG",]$year)){
  coastal_phen[coastal_phen$site_name == "ZACKENBERG" & coastal_phen$year == year_counter ,]$snowmelt_average <- mean_snowmelt[mean_snowmelt$site_name == "ZACKENBERG" & mean_snowmelt$year == year_counter,]$mean_snowmelt
}
# Check whether this has worked
coastal_phen %>% group_by(site_name, year) %>% summarise(snowmelt_mean = first(snowmelt_average)) == mean_snowmelt
# Nice that worked well!


### Calculate variance inflation factors to check whether 
### multicorrelinarity amongst predictors is a problem
# Standart correlation first
cor(cbind(
  coastal_phen$snowmelt,
  coastal_phen$phase_temp,
  coastal_phen$onset_ice_melt
))
# Vifs
vifcor(cbind(
  coastal_phen$snowmelt,
  coastal_phen$phase_temp,
  coastal_phen$onset_ice_melt
))
# Nice no problem! :)
# Snowmelt average
# Standart correlation first
cor(cbind(
  coastal_phen$snowmelt_average,
  coastal_phen$phase_temp,
  coastal_phen$onset_ice_melt
))
# Vifs
vifcor(cbind(
  coastal_phen$snowmelt_average,
  coastal_phen$phase_temp,
  coastal_phen$onset_ice_melt
))

# Standart correlation first
cor(cbind(
  coastal_phen$snowmelt,
  coastal_phen$phase_temp_unmodified,
  coastal_phen$onset_ice_melt
))
# Vifs
vifcor(cbind(
  coastal_phen$snowmelt,
  coastal_phen$phase_temp_unmodified,
  coastal_phen$onset_ice_melt
))
# Nice no problem! :)
# Check correlaiton between temperature variables
cor(cbind(
  coastal_phen$phase_temp,
  coastal_phen$phase_temp_unmodified,
  coastal_phen$spring_temp
))
# Suggests: spring temperatures (May - July) do actually differ quite a bit form the 'responsive' times.

### Within subject mean center predictors
### Subject here is each species - phen phase combination!
# Phase Temperature
coastal_phen$phase_temp_mean <- as.numeric(tapply(coastal_phen$phase_temp, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$phase_temp_rel <- coastal_phen$phase_temp - coastal_phen$phase_temp_mean

# Phase Temperature Unmodified
coastal_phen$phase_temp_unmodified_mean <- as.numeric(tapply(coastal_phen$phase_temp_unmodified, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$phase_temp_unmodified_rel <- coastal_phen$phase_temp_unmodified - coastal_phen$phase_temp_unmodified_mean

# May, June July ("Spring") Temperature
coastal_phen$spring_temp_mean <- as.numeric(tapply(coastal_phen$spring_temp, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$spring_temp_rel <- coastal_phen$spring_temp - coastal_phen$spring_temp_mean


# Snowmelt
coastal_phen$snowmelt_mean <- as.numeric(tapply(coastal_phen$snowmelt, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$snowmelt_rel <- coastal_phen$snowmelt - coastal_phen$snowmelt_mean

# Snowmelt Site Averages
coastal_phen$snowmelt_average_mean <- as.numeric(tapply(coastal_phen$snowmelt_average, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$snowmelt_average_rel <- coastal_phen$snowmelt_average - coastal_phen$snowmelt_average_mean


# Onset melt
coastal_phen$onset_ice_melt_mean <- as.numeric(tapply(coastal_phen$onset_ice_melt, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$onset_ice_melt_rel <- coastal_phen$onset_ice_melt - coastal_phen$onset_ice_melt_mean

# Average spring ice extent
coastal_phen$spring_extent_mean <- as.numeric(tapply(coastal_phen$spring_extent, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$spring_extent_rel <- coastal_phen$spring_extent - coastal_phen$spring_extent_mean


# Year
coastal_phen$year_mean <- as.numeric(tapply(coastal_phen$year, coastal_phen$site_spp_phen, mean)[coastal_phen$site_spp_phen])
coastal_phen$year_rel <- coastal_phen$year - coastal_phen$year_mean
# Year also as factor
coastal_phen$year_fac <- factor(coastal_phen$year)

### Scale variables
coastal_phen$phase_temp_rel_scaled <- scale(coastal_phen$phase_temp_rel, center = F)
coastal_phen$phase_temp_unmodified_rel_scaled <- scale(coastal_phen$phase_temp_unmodified_rel, center = F)
coastal_phen$spring_temp_rel_scaled <- scale(coastal_phen$spring_temp_rel, center = F)
coastal_phen$snowmelt_rel_scaled <- scale(coastal_phen$snowmelt_rel, center = F)
coastal_phen$snowmelt_average_rel_scaled <- scale(coastal_phen$snowmelt_average_rel, center = F)
coastal_phen$onset_ice_melt_rel_scaled <- scale(coastal_phen$onset_ice_melt_rel, center = F)
coastal_phen$spring_extent_rel_scaled <- scale(coastal_phen$spring_extent_rel, center = F)
coastal_phen$year_rel_scaled <- scale(coastal_phen$year_rel, center = F)

# Store atrributes from scaling
scaling_attributes <- data.frame(predictor = c("phase_temp_rel",
                                               "phase_temp_unmodified_rel",
                                               "spring_temp_rel",
                                               "snowmelt_rel",
                                               "snowmelt_average_rel",
                                               "onset_ice_melt_rel",
                                               "onset_spring_extent_rel",
                                               "year_rel"),
                                 scale = c(
                                   attributes(coastal_phen$phase_temp_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$phase_temp_unmodified_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$spring_temp_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$snowmelt_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$snowmelt_average_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$onset_ice_melt_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$spring_extent_rel_scaled)$"scaled:scale",
                                   attributes(coastal_phen$year_rel_scaled)$"scaled:scale"
                                 ))

### Run attribution models ----

## Base model
coastal_phen_base_phase_temp <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                snowmelt + phase_temp + onset_ice_melt + year,
                              random = ~ site_name + site_spp_phen + plot_id +year, 
                              data = coastal_phen, 
                              family = "cengaussian",
                              nitt = 20000)
summary(coastal_phen_base_phase_temp)
save(coastal_phen_base_phase_temp, file = paste0("models/", today, "_coastal_phen_base_phase_temp.Rda"))
# plot(coastal_phen_base)
# Sweet! Hypotheses are met! Let's continue with the more refined models.
# We need to actually account for the whole hierachichal structure.

# First though test base model with phase_temp_unmodified
coastal_phen_base_temp_unmodified <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                snowmelt + phase_temp_unmodified + onset_ice_melt + year,
                              random = ~ site_name + site_spp_phen + plot_id +year, 
                              data = coastal_phen, 
                              family = "cengaussian",
                              nitt = 20000)
summary(coastal_phen_base_temp_unmodified)
save(coastal_phen_base_temp_unmodified, file = paste0("models/", today, "_coastal_phen_base_temp_unmodified.Rda"))

# And spring temperature
coastal_phen_base_unmodified <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                           snowmelt + phase_temp_unmodified + onset_ice_melt + year,
                                         random = ~ site_name + site_spp_phen + plot_id +year, 
                                         data = coastal_phen, 
                                         family = "cengaussian",
                                         nitt = 20000)
summary(coastal_phen_base_unmodified)

## Model with mean centered variables  and random slopes 
# Snomwlet 
coastal_phen_single_snow <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                        snowmelt_rel_scaled + 
                                        year_rel_scaled, 
                                      random = ~ us(1 + snowmelt_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                      data = coastal_phen, 
                                      family = "cengaussian",
                                      nitt = 20000,
                                      pr = T) # Returns full distributions 
summary(coastal_phen_single_snow)
save(coastal_phen_single_snow, file = paste0("models/", today, "_coastal_phen_single_snow.Rda"))

# Snomwlet average
coastal_phen_single_snow_average <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                       snowmelt_average_rel_scaled + 
                                       year_rel_scaled, 
                                     random = ~ us(1 + snowmelt_average_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                     data = coastal_phen, 
                                     family = "cengaussian",
                                     nitt = 20000,
                                     pr = T) # Returns full distributions 
summary(coastal_phen_single_snow_average)
save(coastal_phen_single_snow_average, file = paste0("models/", today, "_coastal_phen_single_snow_average.Rda"))


# Temperature
coastal_phen_single_phase_temp <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                        phase_temp_rel_scaled + 
                                        year_rel_scaled, # add mean in case it explains something
                                      random = ~ us(1 + phase_temp_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                      data = coastal_phen, 
                                      family = "cengaussian",
                                      nitt = 20000,
                                      pr = T) # Returns full distributions 
summary(coastal_phen_single_phase_temp)
save(coastal_phen_single_phase_temp, file = paste0("models/", today, "_coastal_phen_single_phase_temp.Rda"))


# Temperature unmodified
coastal_phen_single_phase_temp_unmodified <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                        phase_temp_unmodified_rel_scaled + 
                                        year_rel_scaled, # add mean in case it explains something
                                      random = ~ us(1 + phase_temp_unmodified_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                      data = coastal_phen, 
                                      family = "cengaussian",
                                      nitt = 20000,
                                      pr = T) # Returns full distributions 
summary(coastal_phen_single_phase_temp_unmodified)
save(coastal_phen_single_phase_temp_unmodified, file = paste0("models/", today, "_coastal_phen_single_phase_temp_unmodified.Rda"))

# Spring Temperature 
coastal_phen_single_spring_temp <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                        spring_temp_rel_scaled + 
                                        year_rel_scaled, # add mean in case it explains something
                                      random = ~ us(1 + spring_temp_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                      data = coastal_phen, 
                                      family = "cengaussian",
                                      nitt = 20000,
                                      pr = T) # Returns full distributions 
summary(coastal_phen_single_spring_temp)
save(coastal_phen_single_spring_temp, file = paste0("models/", today, "_coastal_phen_single_spring_temp.Rda"))

# Sea-ice Onset Melt
coastal_phen_single_onset_melt <- MCMCglmm(cbind(prior_visit, doy) ~ 
                                        onset_ice_melt_rel_scaled + 
                                        year_rel_scaled, # add mean in case it explains something
                                      random = ~ us(1 + onset_ice_melt_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                      data = coastal_phen, 
                                      family = "cengaussian",
                                      nitt = 20000,
                                      pr = T) # Returns full distributions 
summary(coastal_phen_single_onset_melt)
save(coastal_phen_single_onset_melt, file = paste0("models/", today, "_coastal_phen_single_onset_melt.Rda"))

# Sea-ice Spring Extent
coastal_phen_single_spring_extent<- MCMCglmm(cbind(prior_visit, doy) ~ 
                                           spring_extent_rel_scaled + 
                                           year_rel_scaled, # add mean in case it explains something
                                         random = ~ us(1 + spring_extent_rel_scaled + year_rel_scaled):site_spp_phen + site_name + plot_id + year_fac + site_name:year_fac, 
                                         data = coastal_phen, 
                                         family = "cengaussian",
                                         nitt = 20000,
                                         pr = T) # Returns full distributions 
summary(coastal_phen_single_spring_extent)
save(coastal_phen_single_spring_extent, file = paste0("models/", today, "_coastal_phen_single_spring_extent.Rda"))

#plot(coastal_phen_onset_melt_temp)



## Final model call with scaled effects and random slopes for all variables
# Parameter extended prior
a <- 10000
pa_prior <- list(R=list(V=diag(1), nu=0.002),
                 G=list(G1=list(V=diag(5), nu=5, alpha.mu=c(0,0,0,0,0), alpha.V=diag(5)*a),
                        G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                        G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                        G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                        G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

### Full model with all fixed effects, scaled and centered.
# Single run 
coastal_phen_rslopes_full_phen_temp <- MCMCglmm(cbind(prior_visit, doy) ~
                                        snowmelt_rel_scaled +
                                        phase_temp_rel_scaled +
                                        onset_ice_melt_rel_scaled +
                                        year_rel_scaled,
                                      random = ~ us(1 + snowmelt_rel_scaled +
                                                      phase_temp_rel_scaled +
                                                      onset_ice_melt_rel_scaled +
                                                      year_rel_scaled
                                      ):site_spp_phen +
                                        site_name + plot_id + year_fac + site_name:year_fac,
                                      data = coastal_phen,
                                      family = "cengaussian",
                                      nitt = 20000, #1200000,
                                      prior = pa_prior,
                                      pr = T) # Returns full distributions
save(coastal_phen_rslopes_full_phen_temp , file = paste0("models/", today, "_coastal_phen_rslopes_full_phen_temp.Rda"))


# # Run 4 models in parallel to estimate Gelman Rubin criterion
# 
# # Create Cluster
# cluster_one <-  makeCluster(4,type="SOCK", outfile  = paste0(script_path, "clusterMCMCglmm.log"))
# 
# # Export dataset and priors to cluster
# clusterExport(cluster_one, "coastal_phen")
# clusterExport(cluster_one, "pa_prior")
# registerDoSNOW(cluster_one)
# # Run 4 models in parallel
# models <- clusterApply(cluster_one, 1:4, function(i){
#   MCMCglmm::MCMCglmm(cbind(prior_visit, doy) ~
#              snowmelt_rel_scaled +
#              phase_temp_rel_scaled +
#              onset_ice_melt_rel_scaled +
#              year_rel_scaled,
#            random = ~ us(1 + snowmelt_rel_scaled +
#                            phase_temp_rel_scaled +
#                            onset_ice_melt_rel_scaled +
#                            year_rel_scaled
#            ):site_spp_phen +
#              site_name + plot_id + year_fac + site_name:year_fac,
#            data = coastal_phen,
#            family = "cengaussian",
#            nitt = 1200000,
#            prior = pa_prior,
#            pr = T)
# })
# # stop cluster
# stopCluster(cluster_one)
# 
# save(models, file = paste0(script_path, "models/coastal_phen_rslopes_full_4models_2018-03-21.Rda"))
# 
# # Extract SOL (solutions chains) objects from models
# models_sol <- lapply(models, function(m) m$Sol)
# 
# # convert into MCMC.list object
# models_mcmc_list <- as.mcmc.list(models_sol) 
# 
# # prep plot space
# par(mfrow= c(4,2), mar = c(2,2,1,2))
# gelman.plot(models_mcmc_list, auto.layout = F)
# gelman <- gelman.diag(models_mcmc_list)
# save(gelman, file = paste0(script_path, "models/gelmanoutputs_2018-03-21.Rda"))

# Now run for average spring extent
coastal_phen_rslopes_full_spring_extent <- MCMCglmm(cbind(prior_visit, doy) ~
                                                 snowmelt_rel_scaled +
                                                 phase_temp_unmodified_rel_scaled +
                                                 spring_extent_rel_scaled +
                                                 year_rel_scaled,
                                               random = ~ us(1 + snowmelt_rel_scaled +
                                                               phase_temp_unmodified_rel_scaled +
                                                               spring_extent_rel_scaled +
                                                               year_rel_scaled
                                               ):site_spp_phen +
                                                 site_name + plot_id + year_fac + site_name:year_fac,
                                               data = coastal_phen,
                                               family = "cengaussian",
                                               nitt = 20000,#1200000,
                                               prior = pa_prior,
                                               pr = T) # Returns full distributions
save(coastal_phen_rslopes_full_spring_extent, file = paste0("models/", today, "_coastal_phen_rslopes_full_spring_extent.Rda"))

#### Now run with unmodified phase temp
coastal_phen_rslopes_full_phen_temp_unmodified <- MCMCglmm(cbind(prior_visit, doy) ~
                                        snowmelt_rel_scaled +
                                        phase_temp_unmodified_rel_scaled +
                                        onset_ice_melt_rel_scaled +
                                        year_rel_scaled,
                                      random = ~ us(1 + snowmelt_rel_scaled +
                                                      phase_temp_unmodified_rel_scaled +
                                                      onset_ice_melt_rel_scaled +
                                                      year_rel_scaled
                                      ):site_spp_phen +
                                        site_name + plot_id + year_fac + site_name:year_fac,
                                      data = coastal_phen,
                                      family = "cengaussian",
                                      nitt = 20000, #1200000,
                                      prior = pa_prior,
                                      pr = T) # Returns full distributions
save(coastal_phen_rslopes_full_phen_temp_unmodified, file = paste0("models/", today, "_coastal_phen_rslopes_full_phen_temp_unmodified.Rda"))

#### Now run with unmodified spring temp
coastal_phen_rslopes_full_spring_temp <- MCMCglmm(cbind(prior_visit, doy) ~
                                                             snowmelt_rel_scaled +
                                                             spring_temp_rel_scaled +
                                                             onset_ice_melt_rel_scaled +
                                                             year_rel_scaled,
                                                           random = ~ us(1 + snowmelt_rel_scaled +
                                                                           spring_temp_rel_scaled +
                                                                           onset_ice_melt_rel_scaled +
                                                                           year_rel_scaled
                                                           ):site_spp_phen +
                                                             site_name + plot_id + year_fac + site_name:year_fac,
                                                           data = coastal_phen,
                                                           family = "cengaussian",
                                                           nitt = 20000, #1200000,
                                                           prior = pa_prior,
                                                           pr = T) # Returns full distributions
save(coastal_phen_rslopes_full_spring_temp, file = paste0("models/", today, "_coastal_phen_rslopes_full_spring_temp.Rda"))

### Now run with snowmelt site averages
coastal_phen_rslopes_full_snowmelt_averages <- MCMCglmm(cbind(prior_visit, doy) ~
                                                  snowmelt_average_rel_scaled +
                                                  phase_temp_unmodified_rel_scaled +
                                                  onset_ice_melt_rel_scaled +
                                                  year_rel_scaled,
                                                random = ~ us(1 + snowmelt_average_rel_scaled +
                                                                phase_temp_unmodified_rel_scaled +
                                                                onset_ice_melt_rel_scaled +
                                                                year_rel_scaled
                                                ):site_spp_phen +
                                                  site_name + plot_id + year_fac + site_name:year_fac,
                                                data = coastal_phen,
                                                family = "cengaussian",
                                                nitt = 1200000,
                                                prior = pa_prior,
                                                pr = T) # Returns full distributions
save(coastal_phen_rslopes_full_snowmelt_averages  , file = paste0("models/", today, "_coastal_phen_rslopes_full_snowmlet_averages.Rda"))


# Run 4 models in parallel to estimate Gelman Criterion:
# Create Cluster
cluster_one <-  makeCluster(4,type="SOCK", outfile  = "clusterMCMCglmm_unmodified.log")

# Export dataset and priors to cluster
clusterExport(cluster_one, "coastal_phen")
clusterExport(cluster_one, "pa_prior")
registerDoSNOW(cluster_one)
# Run 4 models in parallel
models_unmodified <- clusterApply(cluster_one, 1:4, function(i){
  MCMCglmm::MCMCglmm(cbind(prior_visit, doy) ~
                       snowmelt_rel_scaled +
                       phase_temp_unmodified_rel_scaled +
                       onset_ice_melt_rel_scaled +
                       year_rel_scaled,
                     random = ~ us(1 + snowmelt_rel_scaled +
                                     phase_temp_unmodified_rel_scaled +
                                     onset_ice_melt_rel_scaled +
                                     year_rel_scaled
                     ):site_spp_phen +
                       site_name + plot_id + year_fac + site_name:year_fac,
                     data = coastal_phen,
                     family = "cengaussian",
                     nitt = 1200000,
                     prior = pa_prior,
                     pr = T)
})
# stop cluster
stopCluster(cluster_one)

save(models_unmodified, file = paste0("models/", today, "_coastal_phen_rslopes_unmodified_4models_2018-03-21.Rda"))
# Extract SOL (solutions chains) objects from models
models_sol <- lapply(models_unmodified, function(m) m$Sol)

# convert into MCMC.list object
models_mcmc_list <- as.mcmc.list(models_sol)

# prep plot space
par(mfrow= c(4,2), mar = c(2,2,1,2))
gelman.plot(models_mcmc_list, auto.layout = F)
gelman <- gelman.diag(models_mcmc_list)
save(gelman, file = paste0("models/", today, "_gelman_outputs_unmodified.Rda"))

