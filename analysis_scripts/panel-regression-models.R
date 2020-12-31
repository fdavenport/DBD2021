library(tidyverse)
library(lfe)
library(fixest)
library(stargazer)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

data <- readRDS(state_panel_file) %>%
    ungroup()
data_damages_only <- subset(data, !is.na(damage_value_nozero))

## -----------------------------------------------------------------------------
## MAIN MODELS

linear_model <- felm(damage_value ~ monthly_precip_NORM | STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(linear_model, "../processed_data/regression_models/linear.Rds")

linear_5day <- felm(damage_value ~ precip_5daymax_NORM | STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(linear_5day, "../processed_data/regression_models/linear_5day_model.Rds")

region_model <- felm(damage_value ~ monthly_precip_NORM:region | STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(region_model, "../processed_data/regression_models/region_model.Rds")

region_5day <- felm(damage_value ~ precip_5daymax_NORM:region | STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(region_5day, "../processed_data/regression_models/region_5day_model.Rds")

## -----------------------------------------------------------------------------
## model specification sensitivity

decade_model <- felm(damage_value ~ monthly_precip_NORM:decade | STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(decade_model, "../processed_data/regression_models/decade_model.Rds")


binned_model <- felm(damage_value ~ precip_binstd | STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(binned_model, "../processed_data/regression_models/binned_model.Rds")


quad_model <- felm(damage_value ~ poly(monthly_precip_NORM, 2, raw = TRUE) |
                       STATEYEAR + STATEMONTH,
                   data = data)
saveRDS(quad_model, "../processed_data/regression_models/quad_model.Rds")


prob_model <- felm(damage_binary ~ monthly_precip_NORM | STATEYEAR + STATEMONTH, data = data)
saveRDS(prob_model, "../processed_data/regression_models/prob_model.Rds")


lag_model <- felm(damage_value ~ monthly_precip_NORM:region +
                      monthly_precip_NORM_lag1:region +
                      monthly_precip_NORM_lag2:region | STATEYEAR + STATEMONTH,
                  data = data)
saveRDS(lag_model, "../processed_data/regression_models/lag_model.Rds")


state_model <- felm(damage_value ~ monthly_precip_NORM:STATENAME | STATEYEAR + STATEMONTH,
                    data = data)
saveRDS(state_model, "../processed_data/regression_models/state_model.Rds")


season_model <- felm(damage_value ~  monthly_precip_NORM:region:season |
                         STATEYEAR + STATEMONTH,
                     data = data)
saveRDS(season_model, "../processed_data/regression_models/season_model.Rds")

## -----------------------------------------------------------------------------
## damage transformation sensitivity

main_10mil <- felm(damage_10million ~ monthly_precip_NORM:region | STATEYEAR + STATEMONTH,
                   data = data)
saveRDS(main_10mil, "../processed_data/regression_models/main_10mil_model.Rds")


main_1bil <- felm(damage_1billion ~ monthly_precip_NORM:region | STATEYEAR + STATEMONTH,
                   data = data)
saveRDS(main_1bil, "../processed_data/regression_models/main_1bil_model.Rds")


main_1k <- felm(damage_1k ~ monthly_precip_NORM:region | STATEYEAR + STATEMONTH,
                   data = data)
saveRDS(main_1k, "../processed_data/regression_models/main_1k_model.Rds")


main_asinh <- felm(damage_asinh ~ monthly_precip_NORM:region | STATEYEAR + STATEMONTH,
                   data = data)
saveRDS(main_asinh, "../processed_data/regression_models/main_asinh_model.Rds")


damage_only <- felm(damage_value_nozero ~ monthly_precip_NORM:region | STATEYEAR + STATEMONTH,
                   data = data_damages_only)
saveRDS(damage_only, "../processed_data/regression_models/damage_only_model.Rds")

## -----------------------------------------------------------------------------
## FIXEST MAXIMUM LIKELIHOOD MODELS

poisson_model <- fepois(damage_norm ~ monthly_precip_NORM | STATEYEAR + STATEMONTH,
                        data = data)
saveRDS(poisson_model, "../processed_data/regression_models/poisson_model.Rds")


logit_model <- femlm(damage_binary ~ monthly_precip_NORM | STATEYEAR + STATEMONTH,
             family = "logit",
             data = data)
saveRDS(logit_model, "../processed_data/regression_models/logit_model.Rds")
