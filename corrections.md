## Code error in sea ice data extraction for *Assmann et al. 2019. Global Change Biology* with no consequences for the findings and conclusions 

Jakob J. Assmann j.assmann@bio.au.dk - 14 April 2021

### Summary

- We recently used the sea-ice data extraction code developed for Assmann et al. 2019. Global Change Biology for another project. 
- During the transfer of the code to the new project we spotted a small indexing error in the extractions scripts used in the published manuscript. 
- We resolved the error, re-extracted the data and re-ran the analysis for Assmann et al. 2019. 
- We no noticeable impact on the findings and conclusions presented in the original manuscript. 

---

### Content

1. [The coding error](#the-coding-error)
2. [Little impact on the spring drop in sea-ice predictor](#Little-impact-on-the-spring-drop-in-sea\-ice-predictor)
3. [No noticeable impact on the main findings and conclusions](#No-noticeable-impact-on-the-main-findings-and-conclusions)

---

### The coding error

The coding error occured in the sea-ice data extraction scripts. Specifically, it occured during the step that determines the total amount of daily sea-ice extent in square kilometres within the predefined region for any given site and day. Here, a matrix indexing error caused the size of the individual cells in the NSIDC polar sterographic grid to be read in the wrong way (by columns by row instead of rows by colum):

``` r 
# Original code with error  
for (i in 1:304){
  for (j in 1:448){
    pixel_area[i, j] <- readBin(to_read, integer(), endian = "little", size = 4) / 1000 
  }
}

# Corrected code
for (i in 1:448){
  for (j in 1:304){
    pixel_area_new[i, j] <- readBin(to_read, integer(), endian = "little", size = 4) / 1000 
  }
}
```

The cell sizes were affected as follows:

<img src="data/2020_02_01_quality_control/qc_2021/polar_grid_cell_area_wrong.png" style="zoom:75%;" />

__Figure 1:__ Wrongly loaded cell sizes of the NSIDC polar stereographic grid (North Pole) in square kilometres.



<img src="data/2020_02_01_quality_control/qc_2021/polar_grid_cell_area_correct.png" style="zoom:75%;" />

 __Figure 2:__ Correctly loaded cell sizes of the NSIDC polar stereographic grid  (North Pole) in square kilometres.



The following original scripts were affected:

``` R
data/sea_ice_data/alexfiord/sea_ice_extraction_alexfiord.R # line 65 onwards
data/sea_ice_data/barrow/sea_ice_extraction_barrow.R # line 65 onwards 
data/sea_ice_data/qhi/sea_ice_extraction_qhi.R # line 94 onwards
data/sea_ice_data/zackenberg/sea_ice_extraction_zackenberg.R # line 65 onwards
```

These error was corrected in the following new copies of these scripts:

``` R
data/sea_ice_data/alexfiord/sea_ice_extraction_alexfiord_new.R 
data/sea_ice_data/barrow/sea_ice_extraction_barrow_new.R  
data/sea_ice_data/qhi/sea_ice_extraction_qhi_new.R 
data/sea_ice_data/zackenberg/sea_ice_extraction_zackenberg_new.R
```

---

 ### Little impact on the spring drop in sea-ice predictor

For our anlaysis in Assmann et al. 2019, we used the daily regional sea-ice extent to infer the annual spring drop in sea ice for the four studied regions. The discovered coding error had little effect on the estimted spring drop in regional sea-ice extent.



![](data/2020_02_01_quality_control/qc_2021/diff_sea_ice_drop.png)

__Figure 3__: The coding error had a very small effect on the estimated annual spring drop in sea-ice (Day of Year) for each region. Blue lines and points show the original estimates derieved with the error and red lines and points show the new estimates with the error corrected. 

The spring drop in sea-ice extent can be re-calculated from the corrected data using this script:

``` R
data/dataset_prep_new.R
```

The resulting updated datasets can be found in:

``` R
data/coastal_phen_new.Rda
data/coastal_phen_new.csv
```

---

### No noticeable impact on the main findings and conclusions

The main objective for our study was assess wich of the three predictors (spring snowmelt, spring temperature and spring drop in sea ice extent) were the most meaningful for predicting tundra phenology at four study sites. We re-ran our Baysen hierachical mixed model analysis and found that the corrected onset of spring drop in sea-ice predictor had no noticeable impact on the findings and conclusions presented in the orginal manuscript.



![](analysis/unmodified/coastal_spp_effect_size_plot_centre_legend.png)

__Figure 4:__ **Original figure** showing the estimated effect sizes for the three environmental predictors (snowmelt, temperature and spring drop in sea ice) on the phenology events in Assmann et al. 2019. These results contained the spring drop in sea ice estimates produced with the coding error. 



![](analysis/2021_02_01_new_analysis/coastal_spp_effect_size_plot.png)

__Figure 5:__ Results from the re-analysis using the **corrected** sea ice predictor data show no noticable difference in the estimated effect sizes of the three environmental predictors (snowmelt, temperature and spring drop in sea ice) on the phenology events when compared to the orignal figure from Assmann et al. 2019 (Figure 4).

The updated analysis scripts can be found here:

``` R
analysis/coastal_phen_output_visual_unmodified_new.R # Mixed-model analysis 
analysis/coastal_phen_output_visual_unmodified_new.R # Output visualisation
```

 