library(tidyverse)
library(scales)
library(stargazer)
library(zoo)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

## -----------------------------------------------------------------------------
## READ IN MONTHLY PRECIPITATION DATA
precip_df <- readRDS(monthly_precip_file) %>%
    ## no damages for DC, exclude from analysis
    subset(STATENAME != "District of Columbia")

## calculate mean and sd during baseline period
precip_stats <- precip_df %>%
    subset(year %in% base_period) %>%
    group_by(STATENAME) %>%
    dplyr::summarize(mean_precip = mean(monthly_precip),
              sd_precip = sd(monthly_precip))

precip_df <- left_join(precip_df, precip_stats, by = "STATENAME") %>%
    mutate(monthly_precip_NORM = (monthly_precip - mean_precip)/sd_precip) %>%
    group_by(STATENAME) %>%
    arrange(year, month) %>%
    mutate(monthly_precip_NORM_lag1 = lag(monthly_precip_NORM, 1, order_by = STATENAME),
           monthly_precip_NORM_lag2 = lag(monthly_precip_NORM, 2, order_by = STATENAME)) %>%
    dplyr::select(-mean_precip, -sd_precip)

## save standardized precip data
saveRDS(precip_df, precip_std_file)

## create vector of states to subset damage data
precip_states <- unique(precip_df$STATENAME)

## read daily precip and compute 5-day max for each month
daily_precip_df <- readRDS(daily_precip_file) %>%
    subset(STATENAME != "District of Columbia") %>%
    mutate(year = strftime(as.Date(date), "%Y"),
           month = strftime(as.Date(date), "%m")) %>%
    group_by(STATENAME, year, month) %>%
    arrange(date) %>%
    mutate(precip_5day = rollsum(daily_precip, 5, align = "center", fill = NA)) %>%
    summarize(precip_5daymax = max(precip_5day, na.rm = TRUE))

precip_5day_stats <- daily_precip_df %>%
    subset(year %in% base_period) %>%
    group_by(STATENAME) %>%
    dplyr::summarize(mean_precip_5day = mean(precip_5daymax),
                     sd_precip_5day = sd(precip_5daymax))

daily_precip_df <- daily_precip_df %>%
    left_join(precip_5day_stats, by = c("STATENAME")) %>%
    mutate(precip_5daymax_NORM = (precip_5daymax - mean_precip_5day)/sd_precip_5day) %>%
    dplyr::select(-mean_precip_5day, -sd_precip_5day) %>%
    ungroup() %>%
    mutate(year = as.numeric(year))

saveRDS(daily_precip_df, precip5day_std_file)

## -----------------------------------------------------------------------------
## read in SHELDUS data and aggregate to state-level
damage_data <- read.csv(damage_file) %>%
    dplyr::select(-c(13:21)) %>%  #remove unneeded variables
    dplyr::rename(year = Year) %>%
    mutate(month = str_pad(Month, 2, side = "left", pad = "0")) %>%
    mutate(STATENAME = str_to_title(State.Name)) %>%
    ## only include states in CONUS that are covered by PRISM
    subset(STATENAME %in% precip_states) %>%
    group_by(STATENAME, month, year) %>%
    dplyr::summarize(totaldmg = sum(CropDmg) + sum(PropertyDmg),
              totaldmg_adj2017 = sum(CropDmg.ADJ.2017.) + sum(PropertyDmg.ADJ.2017.))

## -----------------------------------------------------------------------------
## WEALTH DATA

income_data <- read.csv(income_file) %>%
    dplyr::rename(STATENAME = GeoName) %>%
    dplyr::select(-GeoFips, -LineCode) %>%
    mutate(income_100million_current = income_millions_current/100,
           income_10million_current = income_millions_current/10,
           income_billion_current = income_millions_current/1000)

## -----------------------------------------------------------------------------
## COMBINE DATA

data <- left_join(precip_df, daily_precip_df, by = c("STATENAME", "year", "month")) %>%
    left_join(damage_data, by = c("STATENAME", "year", "month")) %>%
    left_join(income_data, by = c("STATENAME", "year")) %>%
    subset(year <= damage_end & year >= damage_start) %>%
    mutate(
        ## NORMALIZED DAMAGE (INCLUDES ZEROS)
        damage_norm = ifelse(totaldmg ==0 | is.na(totaldmg),
                                0,
                             totaldmg/income_100million_current),
        ## LOG NORMALIZED DAMAGE (ZEROS TRANSFORMED)
        damage_value = ifelse(is.na(totaldmg),
                                 log(1),
                       ifelse(totaldmg == 0,
                              log(1),
                              log(totaldmg/income_100million_current))),
        damage_value_nozero = ifelse(damage_value == 0, NA, damage_value),
        damage_10million = ifelse(is.na(totaldmg),
                                     log(1),
                           ifelse(totaldmg == 0,
                                  log(1),
                                  log(totaldmg/income_10million_current))),
        damage_1billion = ifelse(is.na(totaldmg),
                                 log(1),
                          ifelse(totaldmg == 0,
                                 log(1),
                                 log(totaldmg/income_billion_current))),
        damage_1k = ifelse(is.na(totaldmg),
                           log(1000/income_100million_current),
                    ifelse(totaldmg == 0,
                           log(1000/income_100million_current),
                           log(totaldmg/income_100million_current))),
        damage_asinh = ifelse(is.na(totaldmg), asinh(0),
                              asinh(totaldmg/income_100million_current))) %>%
    mutate(region = factor(findNCEIregion(STATENAME), levels = region_order),
           year = as.numeric(year),
           season = factor(ifelse(month %in% c("12", "01", "02"),
                           "DJF",
                    ifelse(month %in% c("03", "04", "05"), "MAM",
                    ifelse(month %in% c("06", "07", "08"), "JJA",
                           "SON"))), levels = season_order)) %>%
    mutate(decade = cut(year, c(1988, 1997, 2007, 2017), include.lowest = TRUE,
                        dig.lab = 4,
                        labels = decade_names)) %>%
    mutate(precip_binstd = cut(monthly_precip_NORM, c(seq(-2.5, 4, 0.5), 7))) %>%
    mutate(damage_binary = ifelse(is.na(totaldmg), 0,
                           ifelse(totaldmg > 0, 1, 0)),
           STATEMONTH = paste0(STATENAME, "_", month),
           STATEYEAR = paste0(STATENAME, "_", year))

saveRDS(data, state_panel_file)
write.csv(data, file = "../processed_data/state_panel_data.csv", row.names = FALSE)

## -----------------------------------------------------------------------------
## NUMBER OF OBSERVED DAMAGE EVENTS

print("Number of state-months with positive flood damage:")
print(nrow(subset(data, totaldmg > 0)))

