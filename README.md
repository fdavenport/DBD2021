Supporting data and code

## Organization of repository: 

* **data**: input data used for analysis (not all raw data is included due to size; see details below)
* **analysis_scripts**: code to read and analyze data
* **processed_data**: output from analysis_scripts (used as input to figure scripts)
* **figure_scripts**: code to create figures 

## data

Datasets not included in the repository are available at the following locations: 
* **PRISM monthly and daily precipitation**: available from the [PRISM Climate
Group](http://www.prism.oregonstate.edu/), Oregon State University
* **CMIP5**: information about CMIP5 and how to access the data can be
  found [here](https://pcmdi.llnl.gov/mips/cmip5/)
* **CLIMDEX HadEX3 historical data**: the CLIMDEX monthly maximum 5-day precipitation ("Rx5day") product is available from [https://www.climdex.org/access/](https://www.climdex.org/access/)
* **CLIMDEX CMIP5 data**: the CLIMDEX CMIP5 datasets are available through [Environment Canada](https://climate-modelling.canada.ca/climatemodeldata/climdex/)
* **Spatial Hazard Events and Losses Database for the U.S.** (SHELDUS): can be accessed with a subscription through the [SHELDUS site](https://cemhs.asu.edu/sheldus)). The aggregated (state-month) flood damages used in the analysis are included in the **processed_data** folder. 
* **NFIP Redacted Claims Dataset**: publicly available from [FEMA](https://www.fema.gov/media-library/assets/documents/180374)

Input data included in this repository: 
* **cb_2017_us_state_20m/**: shapefile of state boundaries from the [Census
  cartographic boundary files](https://www.census.gov/geographies/mapping-files/2017/geo/kml-cartographic-boundary-files.html)
* **income_state_current.csv**: annual state-level income, formatted version of data from the [Bureau of Economic Analysis](https://apps.bea.gov/iTable/index_regional.cfm)
* **fixed_assets_current_dollars.csv**: national net stock of reproducible fixed assets from the [Bureau of Economic Analysis](https://apps.bea.gov/iTable/index_FA.cfm) 
* **combined_home_value_data.csv**: combined data on state-level median home values from the [US Census Housing Tables](https://www.census.gov/hhes/www/housing/census/historic/values.html) for 1990 and 2000, and from the [American Community Survey](https://data.census.gov/cedsci/table?g=0100000US_0400000US29,39,32,31,78,34,33,36,35,38,37,72,30,28,21,20,23,22,66,69,25,24,27,26,60,18,17,19,54,10,53,56,12,55,11,13,16,15,50,51,06,09,08,42,45,01,44,47,02,46,05,49,04,48,41,40&tid=ACSDT1Y2010.B25077&vintage=2010&hidePreview=true&tp=true) for 2010
* **combined_housing_units_data.csv**: combined data on state-level housing units from the US Census for [1990](https://www.census.gov/data/tables/time-series/demo/popest/1990s-housing-units.html), [2000](https://www.census.gov/data/datasets/time-series/demo/popest/intercensal-2000-2010-housing-units.html), and [2010](https://www.census.gov/data/tables/time-series/demo/popest/2010s-total-housing-units.html)
* **WBD_07_HU2_Shape/**: shapefile of the Upper Mississippi River basin from the USGS [Watershed Boundary Dataset](https://www.usgs.gov/core-science-systems/ngp/national-hydrography/watershed-boundary-dataset?qt-science_support_page_related_con=4#qt-science_support_page_related_con)
* **WBD_10_HU2_Shape/**: shapefile of Missouri River basin from the USGS [Watershed Boundary Dataset](https://www.usgs.gov/core-science-systems/ngp/national-hydrography/watershed-boundary-dataset?qt-science_support_page_related_con=4#qt-science_support_page_related_con)
* **billion_dollar_disasters_Flooding.csv**: rainfall and flooding-related events within the [Billion Dollar Weather and Climate Disasters](https://www.ncdc.noaa.gov/billions/overview) 

## analysis_scripts

Some scripts take a while to run - approximate times and memory usage are indicated if run time is >1min.

* **file_paths.R**: file paths for input data and processed data (must be edited for input data not included with repository) 
* **parameters.R**: defines parameter variables used throughout analysis
* **func.R**: defines functions used throughout analysis

1 - Reading PRISM and CMIP precipitation data: 
* **calc-monthly-state-PRISM.R**: calculate monthly timeseries of average precipitation in each state (CALCULATION TIME: ~10min; MEMORY: 16GB)
* **calc-daily-state-PRISM.R**: calculate daily timeseries of average precipitation in each state (CALCULATION TIME: ~1hr 15min; MEMORY: 16GB)
* **calc-gridded-PRISM-timeseries.R**: read regridded (2.5 deg.) monthly PRISM precipitation files and save monthly timeseries (CALCULATION TIME: <5min)
* **read-cmip-2pt5-deg.R**: read 2.5 deg. CMIP5 data and save monthly timeseries (CALCULATION TIME: <5min per model)
* **read-climdex.R**: read 2.5 deg. CLIMDEX data and save monthly timeseries (CALCULATION TIME: ~5min)
* **calc-watersheds-PRISM.R**: calculate monthly timeseries of average precipitation in Missouri River and Upper Mississippi basins (CALCULATION TIME: ~10min)

2 - Combining panel data and running regression models: 
* **combine-state-panel-data.R**: combine monthly flood damages and precipitation data (CALCULATION TIME: <<1min)
* **panel-regression-models.R**: fit regression models (CALCULATION TIME: <<1min)
* **model-bootstrapping.R**: bootstrap regression models (CALCULATION TIME: <2min per model)

3 - Calculating precipitation trends and counterfactual damages:
* **calc-detrended-precip.R**: calculate 10,000 detrended precipitation time series (CALCULATION TIME: each trend time period takes around 15hrs; MEMORY: >12GB)
* **calc-random-precip-trends.R**: calculate random precipitation trends using moving block bootstrap to determine significance of observed trends (CALCULATION TIME: <5min)
* **calc-counterfactual-damages.R**: calculate counterfactual damages for each detrended timeseries and each regression model bootstrap replicate (CALCULATION TIME: each model-time period combination takes about 30min when using 10 parallel tasks; MEMORY: 8GB per task) 

## processed_data

* **monthly_state_precip.Rds**: monthly precipitation time series (state-level)
* **precip_data_std.Rds**: standardized monthly precipitation time series (state-level)
* **daily_state_precip.Rds**: daily precipitation time series (state-level)
* **precip_5daymax.Rds**: standardized monthly maximum 5-day precipitation time series (state-level) 
* **state_panel_data.csv**: combined state-month precipitation and flood damages (csv format)
* **state_panel_data.Rds**: combined state-month precipitation and flood damages (Rds format)
* **watershed_precip_data.Rds**: monthly precipitation time series for Upper Mississippi Basin and Missouri River Basin
* **prism_2pt5deg_precip.Rds**: monthly PRISM precipitation time series on 2.5 degree grid
* **climdex_monthly_Rx5day.Rds**: monthly maximum 5-day precipitation ("Rx5day") time series on 2.5 degree grid
* **2pt5deg_cmip_monthly/**: monthly precipitation time series from cmip models on 2.5 degree grid 
* **2pt5deg_cmip_rx5day/**: monthly max 5-day precipitation ("Rx5day") time series from cmip models on 2.5 degree grid 
* **regression_models/**: regression models (in .Rds format)
* **bootstrapped_models/**: bootstrapped regression model coefficients
* ***detrended_precip/***: detrended state-month precipitation time series (*not included* ~30GB)
* ***counterfactual_damages/***: estimated counterfactual damages (*not included* ~1.5GB) 
* **quantile_trends_1928-2017.Rds**: bootstrapped 50th, 75th, and 95th percentile precipitation trends for each state
* **random_trends.Rds**: random precipitation trends in each state from moving block bootstrap

*files in italics not included in repository due to size*

## figure_scripts
* **theme_func.R**: defines plot themes and plotting functions used in figure scripts
* **fig_1.R**: creates Figure 1
* **fig_2.R**: creates Figure 2
* **fig_3.R**: creates Figure 3 and Figure S6
* **fig_4_5.R**: creates Figure 4, Figure 5, and Figure S7 
* **fig_S1.R**: creates Fig. S1
* **fig_S2.R**: creates Fig. S2
* **fig_S3.R**: creates Fig. S3
* **fig_S4.R**: creates Fig. S4
* **fig_S5.R**: creates Fig. S5
* **fig_S8.R**: creates Fig. S8
* **fig_S9.R**: creates Fig. S9
* **fig_S10.R**: creates Fig. S10
* **fig_S11.R**: creates Fig. S11

## R packages used

* **tidyverse, lfe, fixest, quantreg, raster, ncdf4, sf, velox, RColorBrewer, scales, stringr, ggpubr, stargazer** 

code was written in R 3.5.1 
