# Coastal Phenology Data and Code Repository

## Content
This repository contains the code and data necessary to replicate data analysis, figures, and tables in:

Assmann, Jakob J., Isla H. Myers-Smith, Albert B. Phillimore, Anne D. Bjorkman, Richard E. Ennos, Janet S. Prevéy, Greg H.R. Henry, Niels M. Schmidt and Robert D. Hollister. **In press**. ***Local snowmelt and temperature – but not regional sea-ice – explain variation in spring phenology in coastal Arctic tundra***. *Global Change Biology*.

## Contact
Jakob Assmann 

Email: jakobjassmann [at] gmail . com

Website: [jakobjassmann.wordpress.com](https://jakobjassmann.wordpress.com/)
# Data preparation

All data for the phenological observations and environmental predictors is already publicly available. The only exception is the Zackenberg plot-level phenology observations, which we have included with the data preparation scripts in this repostiorty. 

The data preparation, cleaning and assembly scripts can be found in the following locations:
```
/data/phenology_data
/data/temperature_data
/data/sea_ice_data
```
Links to the raw data sets are provided in the readme files within each of the folder. Please follow these links for the respective data sharing and usage agreements of the respective datasets. 

The fowllowing script compiles the data prepared with the scripts above into a single dataset for later use in the analysis (see next section):
```
/data/dataset_prep.R
```

# Data 
Interval censored phenology observation per site and plot with associated environmental predictor variables can be found in the following location:
```
/data/coastal_phen.Rda
/data/coastal_phen.csv
```
The content of the Rda and CSV files is identical.

# Analysis scripts

```
```