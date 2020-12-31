library(tidyverse)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

## BATCH PARAMETERS
args <- commandArgs(trailingOnly = TRUE)
ARRAY_ID <- as.numeric(args[1]) ## job number within array - used to determine which sensitivity combination to calculate
N <- as.numeric(args[2])
modelcode <- args[3]
trend_range <- args[4]
LAG <- grepl("lag", modelcode)
QUAD <- grepl("quad", modelcode)
if(QUAD){
    POLY <- TRUE
} else {
    POLY <- FALSE
}
## -----------------------------------------------------------------------------
## INPUT DATA
damage_data <- readRDS(state_panel_file) %>%
    dplyr::select(STATENAME, year, month, region,
                  season, totaldmg_adj2017) %>%
    ## only calculate counterfactual damages for months with observed damage
    subset(totaldmg_adj2017 > 0)

## -----------------------------------------------------------------------------
## read in detrended precip
detrended_precip_list <- lapply(paste0("../processed_data/detrended_precip/detrended_precip_",
                                         trend_range, "_",
                                         ((ARRAY_ID-1)*N+1):((ARRAY_ID)*N), ".Rds"),
                                readRDS)

## -----------------------------------------------------------------------------
## read and format in model coefficients
model_coeff_boot <- read_boot_results(paste0("../processed_data/bootstrapped_models/", modelcode,
                                            "_bootstrap_1.Rds"),
                                     n_var = 3,
                                     var_codes = c("precip", "region", "season"), POLY) %>%
    dplyr::rename(coeff = Estimate) %>%
    dplyr::select(-`Std. Error`)

## remove season column if all NAs
if(sum(!is.na(model_coeff_boot$season)) == 0){
    model_coeff_boot <- dplyr::select(model_coeff_boot, -season)
}
if(sum(!is.na(model_coeff_boot$region)) == 0){
    model_coeff_boot <- dplyr::select(model_coeff_boot, -region) %>%
        crossing(region = unique(damage_data$region))
}
if(POLY){
    model_coeff_boot <- pivot_wider(model_coeff_boot, names_from = poly_order,
                               names_prefix = "poly", values_from = coeff)
}

subresult <- vector('list', length = length(detrended_precip_list))

for(i in seq_along(detrended_precip_list)){
    if(i %% 10 == 0) print(i)

    if(LAG){
        monthly_damages <- detrended_precip_list[[i]] %>%
            dplyr::select(STATENAME, year, month, precip_diff) %>%
            ## CALCULATE LAGGED PRECIP FOR EACH STATE
            group_by(STATENAME) %>%
            arrange(year, month) %>%
            mutate(monthly_precip_NORM_lag1_diff = lag(precip_diff, 1,
                                                       order_by = STATENAME),
                   monthly_precip_NORM_lag2_diff = lag(precip_diff, 2,
                                                       order_by = STATENAME)) %>%
            ungroup() %>%
            right_join(damage_data, by = c("STATENAME", "year", "month")) %>%
            ## SPREAD PRECIP VARIABLES AND CALCULATE CHANGE FACTOR
            dplyr::rename(monthly_precip_NORM_diff = precip_diff) %>%
            gather(precip, detrended_diff, contains("_diff")) %>%
            mutate(precip = gsub("_diff", "", precip)) %>%
            left_join(model_coeff_boot) %>%
            mutate(beta_prod = coeff*(-detrended_diff)) %>%
            dplyr::select(-detrended_diff, -coeff) %>%
            pivot_wider(names_from = precip, values_from = beta_prod) %>%
            mutate(change_factor = exp(monthly_precip_NORM +
                                       monthly_precip_NORM_lag1 +
                                       monthly_precip_NORM_lag2),
                   cf_dmg = totaldmg_adj2017*change_factor)
    } else if (QUAD){
        monthly_damages <- detrended_precip_list[[i]] %>%
            dplyr::select(STATENAME, year, month, monthly_precip_NORM, precip_diff) %>%
            right_join(damage_data, by = c("STATENAME", "year", "month")) %>%
            left_join(model_coeff_boot) %>%
            ## CALCUALTE BETA FROM QUADRATIC
            mutate(beta = poly1 + poly2*monthly_precip_NORM,
                   change_factor = exp(beta*(-precip_diff)),
                   cf_dmg = totaldmg_adj2017*change_factor)
    } else {
        monthly_damages <- detrended_precip_list[[i]] %>%
            dplyr::select(STATENAME, year, month, precip_diff) %>%
            right_join(damage_data, by = c("STATENAME", "year", "month")) %>%
            left_join(model_coeff_boot) %>%
            mutate(beta_prod = coeff*(-precip_diff),
                   change_factor = exp(beta_prod),
                   cf_dmg = totaldmg_adj2017*change_factor)

    }

    cmltive_dmgs <- monthly_damages %>%
        group_by(boot_id) %>%
        summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
                  cf_dmg= sum(cf_dmg, na.rm = TRUE))

    subresult[[i]] <- cmltive_dmgs
}

result <- bind_rows(subresult) %>%
    mutate(model_code = modelcode,
           trend_range = trend_range)

saveRDS(result, paste0("../processed_data/counterfactual_damages/damage_sensitivity_",
                       trend_range, "_", modelcode, "_", ARRAY_ID, ".Rds"))

