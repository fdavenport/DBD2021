## -----------------------------------------------------------------------------
## PARAMETERS
## -----------------------------------------------------------------------------

## reading prism precipitation data:
prism_monthly_startyear <- 1895
prism_monthly_endyear <- 2017
prism_daily_years <- 1981:2017

# time period for analysis of damages:
damage_start <- 1988
damage_end <- 2017

# time period for precipitation trend calculation:
precip_start <- 1928
precip_end <- 2017
precip_range <- paste0(precip_start, "-", precip_end)

base_period <- 1986:2005 ## time period for standardizing precipitation

## -----------------------------------------------------------------------------
## QUANTILE TREND CALCULATIONS
quants <- c(0.5, 0.75, 0.95)
boot_block_size <- 12 ## months
rq_nboot <- 1000 ## number of bootstraps for random trend calculation

## -----------------------------------------------------------------------------
## COUNTERFACTUAL TIME SERIES CALCULATION
q_breaks <- seq(.1, .9, .1)
q_mid <- seq(0.05, .95, .1)
d_prop <- .2 ## 'delete-d proportion': proportion of years to delete for trend uncertainty

## -----------------------------------------------------------------------------
## COUNTERFACTUAL SENSITIVITY TESTS

## precip. trend sensitivity:
ppt_trend_starts <- c(1898, 1908, 1918, 1928, 1938, 1948, 1958)
ppt_trend_end <- 2017
ppt_trend_ranges <- paste0(ppt_trend_starts, "-", ppt_trend_end)

## model specification sensitivity (named vector of model codes):
model_sensitivity_list <- c("regional\nmodel (main)" = "region_model",
                            "+ season\ninteraction" = "season_model",
                            "+ lagged\nprecipitation" = "lag_model",
                            "linear\n(all regions)" = "linear",
                            "quadratic\n(all regions)" = "quad_model")
model_sensitivity_labels <- names(model_sensitivity_list)

## dependent variable transformation sensitivity (named vector of model codes):
transform_sensitivity_list <- c("$1/10M\nincome" = "main_10mil_model",
                            "$1/100M\nincome\n (main)" = "region_model",
                            "$1/$1B\nincome" = "main_1bil_model",
                            "$1000" = "main_1k_model",
                            "excluded" = "damage_only_model",
                            "inverse\nhyperbolic\nsine\ntransformation" = "main_asinh_model")
transform_sensitivity_labels <- names(transform_sensitivity_list)

## -----------------------------------------------------------------------------
## CMIP5 ANALYSIS:
## ipcc simulations (for subsetting analysis)
ipcc_runs <- c("ACCESS1-0_r1i1p1", "ACCESS1-3_r1i1p1",
               "bcc-csm1-1_r1i1p1", "BNU-ESM_r1i1p1","CanESM2_r1i1p1",
               "CCSM4_r1i1p1", "CESM1-BGC_r1i1p1",
               "CESM1-CAM5_r1i1p1", "CMCC-CM_r1i1p1",
               "CMCC-CMS_r1i1p1", "CNRM-CM5_r1i1p1",
               "CSIRO-Mk3-6-0_r1i1p1", "EC-EARTH_r8i1p1",
               "FGOALS-g2_r1i1p1", "FIO-ESM_r1i1p1",
               "GFDL-CM3_r1i1p1", "GFDL-ESM2G_r1i1p1",
               "GFDL-ESM2M_r1i1p1", "GISS-E2-H_r1i1p1",
               "GISS-E2-R_r1i1p1", "HadGEM2-AO_r1i1p1",
               "HadGEM2-CC_r1i1p1", "HadGEM2-ES_r2i1p1",
               "inmcm4_r1i1p1", "IPSL-CM5A-LR_r1i1p1",
               "IPSL-CM5A-MR_r1i1p1", "IPSL-CM5B-LR_r1i1p1",
               "MIROC5_r1i1p1", "MIROC-ESM_r1i1p1",
               "MIROC-ESM-CHEM_r1i1p1", "MPI-ESM-LR_r1i1p1",
               "MPI-ESM-MR_r1i1p1", "MRI-CGCM3_r1i1p1",
               "NorESM1-M_r1i1p1", "NorESM1-ME_r1i1p1")

## -----------------------------------------------------------------------------
## FACTOR VARIABLE LEVELS AND LABELS

## REGIONS
region_order <- c("Northwest", "West", "Southwest", "Northern Rockies Plains",
                  "South", "Upper Midwest", "Central", "Southeast", "Northeast",
                  "otherNCEI")
## STATES
state_order <- data.frame(STATENAME = state.name) %>%
    mutate(abb = state.abb[match(STATENAME, state.name)],
           region = factor(findNCEIregion(STATENAME), levels = region_order)) %>%
    arrange(region) %>%
    ungroup()
state_order <- state_order %>%
    mutate(abb = factor(abb, levels = state_order$abb))

## SEASONS
season_order <- c("DJF", "MAM", "JJA", "SON")

## DECADES
dec1 <- seq(1988, 2008, 10)
dec2 <- seq(1997, 2017, 10)
decade_names <- paste0(dec1, "-\n", dec2)




