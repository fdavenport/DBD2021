library(tidyverse)
library(quantreg)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

## -----------------------------------------------------------------------------
## calculate precip trends using delete-D method
## all states recalculated together (states are not independent)

args <- commandArgs(trailingOnly = TRUE)
ARRAY_ID <- as.numeric(args[1]) ## job number within array - used to set seed and in outfile name

## INPUT DATA:
precip_df <- readRDS(precip_std_file)
## start one year earlier than damage period to calculate lag months later on
precip_df_damage_period <- subset(precip_df, year %in% (damage_start-1):damage_end)

## -----------------------------------------------------------------------------
states <- unique(precip_df$STATENAME)
nboot <- 10000

set.seed(123*ARRAY_ID)
k <- ARRAY_ID

subresult <- vector('list', length = nboot)
quant_trends <- vector('list', length = nboot)

print("Trend starting year:")
print(ppt_trend_starts[k])
for(i in 1:nboot){
    if(i %% 100 == 0) print(i)
    ## sample years to include in trend calculations
    include_years <- sample(ppt_trend_starts[k]:ppt_trend_end,
                            length(ppt_trend_starts[k]:ppt_trend_end)*(1-d_prop),
                            replace = FALSE)
    precip_subset <- subset(precip_df, year %in% include_years)

    ## calculate quantile breaks for years within damage period
    quantile_breaks <- calc_bin_breaks(precip_subset, q_breaks,
                                       allyears = (damage_start-1):damage_end)

    ## calculate trends for each quantile bin
    quantile_trends <- calc_bin_trends(precip_subset, q_mid) %>%
        mutate(trend_start = ppt_trend_starts[k])

    ## calculate detrended monthly precip within damage period
    subresult[[i]] <- detrend_precip(precip_df_damage_period,
                                     quantile_breaks,
                                     quantile_trends,
                                     trend_start = ppt_trend_starts[k]) %>%
        mutate(boot_id = i)

    saveRDS(subresult[[i]],
            file = paste0("../processed_data/detrended_precip/detrended_precip_",
                          ppt_trend_ranges[k], "_", i, ".Rds"))

    ## calculate trends at quant levels for extended data figure
    quant_trends[[i]] <- calc_bin_trends(precip_subset, quants) %>%
        mutate(trend_start = ppt_trend_starts[k],
               boot_id = i)

    if(i %% 250 == 0) {
        saveRDS(bind_rows(quant_trends),
                paste0("../processed_data/quantile_trends_", ppt_trend_ranges[k], ".Rds"))
    }

}

saveRDS(bind_rows(quant_trends),
        paste0("../processed_data/quantile_trends_", ppt_trend_ranges[k], ".Rds"))
