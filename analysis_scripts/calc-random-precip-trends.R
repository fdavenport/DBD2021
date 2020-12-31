library(tidyverse)
library(quantreg)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

precip_df <- readRDS(precip_std_file)

##-----------------------------------------------------------------------------
## bootstrap precip trends using moving block bootstrap to estimate significance
## assuming states are independent
states <- unique(precip_df$STATENAME)

set.seed(123)

state_results <- vector('list', length(states))

for(k in seq_along(states)){
    print(k)
    print(states[k])
    tmpdat <- subset(precip_df, STATENAME == states[k] & year >= precip_start &
                                year <= precip_end)

    state_results[[k]] <- as.data.frame(mbb_rq_ts(tmpdat$year, tmpdat$monthly_precip_NORM,
                                                  tau = quants,
                                                  block_size = boot_block_size,
                                                  nboot = rq_nboot)) %>%
        mutate(STATENAME = states[k],
               range = precip_range)
}

random_trends <- bind_rows(state_results)

saveRDS(random_trends, random_trend_file)

