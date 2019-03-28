# Coastal Phenology Data and Code Repository

## Content
This repository contains the code and data necessary to replicate data analysis, figures, and tables in:

Assmann, Jakob J., Isla H. Myers-Smith, Albert B. Phillimore, Anne D. Bjorkman, Richard E. Ennos, Janet S. Prevéy, Greg H.R. Henry, Niels M. Schmidt and Robert D. Hollister. In press. ***Local snowmelt and temperature – but not regional sea-ice – explain variation in spring phenology in coastal Arctic tundra***. Global Change Biology.

## Contact
Jakob J. Assmann 

Email: j.assmann [at] ed.ac.uk

Website: [jakobjassmann.wordpress.com](https://jakobjassmann.wordpress.com/)

## Data & code usage guidelines and license 
### Data
All data for the phenological observations and environmental predictors is already publicly available. Links to the datasets are listed below. Please refer to the relevant data usage guidleines of each datasets.

Phenological observations:
- [Phenological Observations and Snowmlet (Alexandra Fiord and Utqiaġvik - Barrow)](www.polardata.ca/pdcsearch/PDCSearchDOI.jsp?doi_id=12722) (Prevey et al. 2017)
- Additonal observations for Qikqitaruk: [Phenological Observations and Snowmelt (Qikiqtaruk - Herschel Island)](https://github.com/ShrubHub/QikiqtarukHub/blob/master/data/qhi_phen_with_before_2017.csv) (Myers-Smith et al. 2019)
- Phenological observations for Zackenberg were provided by the Greenland Ecosystem Monitoring Programme. Data available at: http://data.g-e-m.dk A pre-formatted version of this data is included in the above PDC dataset and this repository includes a version of these data with additional plot-level observations.

Environmental predictors:
- Snowmelt observations are available via the above phenology observation datasets.
- Temperature observations for Alexandra Fiord are included in this repository. These data were provided by Anne Bjorkman and Greg Henry (Bjorkman et al., 2015). Please contact authors for guidance on data usage.
- Temperature observations for Utqiaġvik - Barrow from the NOAA Earth System Research Laboratory Utqiaġvik Global Monitoring Division. Data available at: https://www.esrl.noaa.gov/gmd/obop/brw/ (NOAA ESRL Global Monitoring Division, 2018)
- Temperature observations for Qiqikqtaruk from Environment Canada Qikiqtaruk - Herschel Island weather station (ID 1560) Daily gap-filled with Environment Canada Komakuk (ID 10822) - see methods. Data available at:
http://climate.weather.gc.ca/historical_data/search_historic_data_e.html
- Temperature observations for Zackenberg provided by the Greenland Ecosystem Monitoring Programme. Data available at: http://data.g‐e‐m.dk.
- Sea-Ice observations were obtained from the [NOAA/NSIDC Climate Data Record v3 Passive Microwave Sea Ice Contentrations](https://nsidc.org/data/G02202) (Meier et al., 2017; Peng, Meier, Scott, & Savoie, 2013)

### Code
All code provided for data preparation and analysis is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/). In accordance with the license the code is available to be shared and adapted, but requires attribution to the authors, e.g. through citation of the above manuscript, and indications where changes were made. 
Although not mandatory, we additionally suggest that code users contact and collaborate with contributors should the code form a substantial proportion of a particular publication or analysis.

# Data preparation
The data preparation, cleaning and assembly scripts can be found in the following locations:
```
/data/phenology_data
/data/temperature_data
/data/sea\_ice\_data
```
The following script compiles the data prepared with the scripts above into a single dataset for later use in the analysis (this script also produces Figure S4):
```
/data/dataset_prep.R
```

# Data 
Interval censored phenology observation per site, species and plot with associated environmental predictor variables can be found in the following location:
```
/data/coastal_phen.Rda
/data/coastal_phen.csv
```
The content of the .Rda and .csv files is identical.

*Important:* Please note that we provide this summarised data for archival purposes only. If you intent to use the phenological observations in this dataset please see data usage and guiance for the raw data sets described above. 

# Analysis scripts
The following scripts can be used to conduct the analysis and produce the figures:
```
Map figure (Figure 1):
/analysis/coastal\_site\_map.r

Phenology trends (statistical models, Figure 2 and S7,  and Table S6):
/analysis/coastal_phentrends.r

Environmental predictor trends (statistical models, Figure 3, Table S8):
/analysis/coastal\_phen\_env\_trends\_unmodified.r

Prediction analysis (statistical models, Table S10):
/analysis/coastal\_phen\_attribution.r

Prediction analysis (Figures 4, S9, Table S3, S11):
/analysis/coastal_phen_output_visual_unmodified.r		
(Using average snowmelt at a site instead of plot-level snomwlet)
/analysis/coastal_phen_output_visual_snowmelt_avg.r	
```
*Please note:* Due to a developmental legacy the data set and analysis contain two versions of the temperature variable:
- Version 1: 'temperature_unmodified' Mean of the daily spring temperature in the period from first snowmelt on record to 75% of the phenology events occuring for a given site-species-phenological event combination. *This is the temperature variable used in the analysis presented in the paper.*
- Version 2: 'temperature' Mean of the daily spring temperature in the period from two weeks prior first snowmelt on record to 75% of the phenology events occuring for a given site-species-phenological event combination. This temperature variable was not used in the analysis presented in the paper.

# References
Bjorkman, A. D., Elmendorf, S. C., Beamish, A. L., Vellend, M., & Henry, G. H. R. (2015). Contrasting effects of warming and increased snowfall on Arctic tundra plant phenology over the past two decades. Global Change Biology, 21, 4651–4661. https://doi.org/10.1111/gcb.13051

NOAA ESRL Global Monitoring Division. (2018). Meteorology measurements at Barrow Atmospheric Baseline Observatory,  Alaska. Boulder, Colorado, USA: Compiled by the Observatory Operations Group. National Oceanic and Atmospheric Administration (NOAA), Earth System Research Laboratory  (ESRL), Global Monitoring Division (GMD).

Meier, W. N., Fetterer, F., Savoie, M., Mallory, S., Duerr, R., & Stroeve, J. (2017). NOAA/NSIDC Climate Data Record of Passive 
Microwave Sea Ice Concentration, Version 3. NOAA/NSIDC daily sea Ice CDR. Retrieved from NSIDC: National Snow and Ice Data Center website: https://doi.org/10.7265/N59P2ZTG

Myers-Smith, I. H., Grabowski, M. M., Thomas, H. J. D., Angers-Blondin, S., Daskalova, G. N., Bjorkman, A. D., … Eckert, C. D. (2019). Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change. Ecological Monographs, e01351. https://doi.org/10.1002/ecm.1351

Peng, G., Meier, W. N., Scott, D. J., & Savoie, M. H. (2013). A long-term and reproducible passive microwave sea ice concentration data record for climate studies and monitoring. Earth System Science Data, 5(2), 311–318. https://doi.org/10.5194/essd-5-311-2013

Prevéy, J. S., Vellend, M., Rüger, N., Hollister, R. D., Bjorkman, A. D., Myers-Smith, I. H., … Rixen, C. (2017). Greater temperature sensitivity of plant phenology at colder sites: implications for convergence across northern latitudes. Global Change Biology, 23(7), 2660–2671. https://doi.org/10.1111/gcb.13619

