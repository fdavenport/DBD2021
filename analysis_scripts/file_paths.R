## blank paths should be edited based on location of input data

## -----------------------------------------------------------------------------
## PRISM DATA

prism_monthly_path <- ""
prism_daily_path <- ""

## monthly precipitation data on 2.5 x 2.5 degree grid:
prism_2pt5deg_path <- ""

## -----------------------------------------------------------------------------
## CMIP5 MONTHLY DATA
## on 2.5 x 2.5 degree grid (each forcing in a separate directory):
cmip_hist_path <- ""
cmip_rcp85_path <- ""
cmip_rcp26_path <- ""

## -----------------------------------------------------------------------------
## CLIMDEX Rx5day DATA (on 2.5 degree grid):

## HadEX3 monthly Rx5day data:
climdex_obs_rx5day_file <- ""
## CMIP5 monthly Rx5day data:
climdex_cmip_path <- ""

## -----------------------------------------------------------------------------
## SHELDUS DATA (damages for events with flooding):
damage_file <- ""

## NFIP CLAIMS DATA:
nfip_claims_file <- ""

## BILLION DOLLAR DISASTER DATA:
bdd_file <- "../data/billion_dollar_disasters_Flooding.csv"

## -----------------------------------------------------------------------------
## GEOGRAPHIC BOUNDARY FILES:
state_boundary_file <- "../data/cb_2017_us_state_20m/cb_2017_us_state_20m.shp"
mo_basin_file <- "../data/WBD_10_HU2_Shape/WBDHU2.shp"
um_basin_file <- "../data/WBD_07_HU2_Shape/WBDHU2.shp"

## -----------------------------------------------------------------------------
## WEALTH/EXPOSURE DATA
assets_file <- "../data/fixed_assets_current_dollars.csv"
income_file <- "../data/income_state_current.csv"
home_value_state_file <- "../data/combined_home_value_data.csv"
housing_units_file <-"../data/combined_housing_unit_data.csv"

## -----------------------------------------------------------------------------
## PROCESSED DATA FILE PATHS:

## monthly state precipitation:
monthly_precip_file <- "../processed_data/monthly_state_precip.Rds"

## standardized monthly state precipitation:
precip_std_file <- "../processed_data/precip_data_std.Rds"

## daily state precipitation
daily_precip_file <- "../processed_data/daily_state_precip.Rds"

## standardized max 5-day state precipitation:
precip5day_std_file <- "../processed_data/precip_5daymax.Rds"

## state panel data:
state_panel_file <- "../processed_data/state_panel_data.Rds"

## random state precipitation trends (using moving block bootstrap):
random_trend_file <- "../processed_data/random_trends.Rds"

## monthly prism timeseries (2pt5 degree grid):
monthly_prism_2pt5_file <- "../processed_data/prism_2pt5deg_precip.Rds"

## rx5day time series (observed, 2pt5 degree grid):
rx5day_obs_file <- "../processed_data/climdex_monthly_Rx5day.Rds"

## cmip monthly precip timeseries, 2pt5 degree grid (directory):
monthly_cmip_path <- "../processed_data/2pt5deg_cmip_monthly"

## cmip Rx5day timeseries, 2pt5 degree grid (directory):
rx5day_cmip_path <- "../processed_data/2pt5deg_cmip_rx5day"

## Missouri and Upper Mississippi watershed precip timeseries:
watershed_precip_file <- "../processed_data/watershed_precip_data.Rds"
