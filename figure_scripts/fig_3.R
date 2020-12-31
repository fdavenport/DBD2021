library(tidyverse)
library(quantreg)
library(ggpubr)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
damage_data <- readRDS(state_panel_file)
precip_data <- readRDS(precip_std_file)

## -----------------------------------------------------------------------------
## read in regional coefficient data
region_coeff <- format_coeff(readRDS("../processed_data/regression_models/region_model.Rds"),
                              n_var = 2,
                             var_codes = c("precip", "region")) %>%
    dplyr::select(region, precip, coeff = Estimate)

## -----------------------------------------------------------------------------
## read in bootstrapped damage results
damage_boot_df <- bind_rows(lapply(
    list.files(path = "../processed_data/counterfactual_damages/",
               pattern = paste0("damage_sensitivity_1928-2017_region_model_"),
               full.names = TRUE), readRDS))

damage_diff_limits <- damage_boot_df %>%
    summarize(obs_damage = unique(obs_damage),
              low = obs_damage - quantile(cf_dmg, 0.975),
              high = obs_damage - quantile(cf_dmg, 0.025)) %>%
    mutate(year = 2017)

## -----------------------------------------------------------------------------
## CALCULATE DETRENDED PRECIP AND DAMAGES

quantile_breaks <- calc_bin_breaks(subset(precip_data,
                                          year >= precip_start &
                                          year <= precip_end),
                                   q_breaks,
                                   allyears = precip_start:precip_end)
quantile_trends <- calc_bin_trends(subset(precip_data,
                                          year >= precip_start &
                                          year <= precip_end), q_mid)
detrended_precip <- detrend_precip(subset(precip_data,
                                          year >= precip_start &
                                          year <= precip_end), quantile_breaks,
                                   quantile_trends,
                                   trend_start = precip_start)

cumulative_damages <- calc_cf_damage(detrended_precip, damage_data, region_coeff) %>%
    group_by(year) %>%
    summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
              cf_dmg = sum(cf_dmg, na.rm = TRUE)) %>%
    arrange(year) %>%
    mutate(total_obs = cumsum(obs_damage),
           cumulative_cf = cumsum(cf_dmg),
           total_diff_mid = total_obs - cumulative_cf)

damage_text <- data.frame(year = 2017.8,
                 value = as.numeric(c(cumulative_damages[nrow(cumulative_damages),"total_obs"],
                           cumulative_damages[nrow(cumulative_damages),"total_diff_mid"],
                           damage_diff_limits[1,"low"],
                           damage_diff_limits[1,"high"]))) %>%
    mutate(label = paste0("$", round(value/1e9,0), "B"))

print("2017 cumulative damages:")
print(subset(cumulative_damages, year == 2017))
print("2017 cumulative percent:")
print(paste0(subset(cumulative_damages, year == 2017)$total_diff_mid /                                            subset(cumulative_damages, year == 2017)$total_obs*100, "%"))
print("2017 damage bounds:")
print(damage_diff_limits)
print("2017 percent bounds:")
print(paste0(damage_diff_limits$low/damage_diff_limits$obs_damage*100, "%, ",
             damage_diff_limits$high/damage_diff_limits$obs_damage*100, "%"))

## -----------------------------------------------------------------------------
## FIGURE 3A:

damage_col <- c("1obs" = "slategray4", "2cf" = "forestgreen")

fig_3A <- ggplot(cumulative_damages[-1,],
                               aes(x = year, y = total_diff_mid/1e9)) +
    geom_line(aes(col = "1obs", y = total_obs/1e9)) +
    geom_line(aes(col = "2cf")) +
    geom_area(fill = "forestgreen", alpha = 0.3) +
    geom_ribbon(fill = "slategray4", alpha = 0.2, aes(ymin = total_diff_mid/1e9, ymax = total_obs/1e9)) +
    geom_point(col = "slategray4", data = subset(cumulative_damages, year == 2017),
               aes(y = total_obs/1e9)) +
    geom_point(col = "forestgreen", data = subset(cumulative_damages, year == 2017)) +
    geom_errorbar(data = damage_diff_limits, col = "forestgreen",
                  aes(ymin = low/1e9, ymax = high/1e9,
                      y = NULL)) +
    geom_text(data = damage_text, size = 2, hjust = 0,
              aes(x = year, y = value/1e9, label = label)) +
    scale_color_manual(values = damage_col,
                       labels = c("Cumulative historical flood damages",
                                  "Estimated portion due to \nprecipitation change")) +
    scale_y_continuous(name = "Billions (2017$)", expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0), limits = c(1988, 2021),
                       breaks = seq(1990, 2010, 10)) +
    ggtitle("Contribution of historical precipitation change to US flood damages") +
    theme_pub +
    theme(axis.title.x = element_blank(),
          legend.position = c(.38, .8),
          legend.title = element_blank(),
          legend.background = element_blank(),
          plot.title = element_text(hjust = 1, face = "bold")) +
    coord_cartesian(ylim = c(0,220))


## -----------------------------------------------------------------------------
## SENSITIVITY TO MODEL SPECIFICATION
## -----------------------------------------------------------------------------

base_precip <- subset(precip_data,
                      year >= precip_start & year <= precip_end)
quantile_breaks <- calc_bin_breaks(base_precip, q_breaks,
                                   allyears = (damage_start-1):damage_end)
quantile_trends <- calc_bin_trends(base_precip, q_mid)
base_precip_detrended <- detrend_precip(base_precip,
                                       quantile_breaks,
                                       quantile_trends,
                                       trend_start = precip_start) %>%
    right_join(damage_data)

model_sensitivity_results <- lapply(
    model_sensitivity_list,
    function(x)
        bind_rows(lapply(list.files(path = "../processed_data/counterfactual_damages/",
                          pattern = paste0("damage_sensitivity_1928-2017_", x),
                          full.names = TRUE), readRDS)))


diff_limits_model <- vector('list', length = length(model_sensitivity_results))

for(k in seq_along(model_sensitivity_list)){
    LAG <- grepl("lag", model_sensitivity_list[k])
    QUAD <- grepl("quad", model_sensitivity_list[k])
    POLY <- ifelse(QUAD, TRUE, FALSE)
    ## -------------------------------------------------------------------------
    ## calculate midpoint
    model_coeff <- format_coeff(readRDS(paste0("../processed_data/regression_models/",
                                       model_sensitivity_list[k],
                                       ".Rds")),
                                n_var = 3,
                                var_codes = c("precip", "region", "season"), POLY) %>%
        dplyr::rename(coeff = Estimate) %>%
        dplyr::select(-`Std. Error`)

     ## remove season column if all NAs
    if(sum(!is.na(model_coeff$season)) == 0){
        model_coeff <- dplyr::select(model_coeff, -season)
    }
    if(sum(!is.na(model_coeff$region)) == 0){
        model_coeff<- dplyr::select(model_coeff, -region) %>%
            crossing(region = unique(damage_data$region))
    }

    damage_mid <- calc_cf_damage(base_precip_detrended, damage_data,
                                 model_coeff, LAG, QUAD) %>%
        summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
                  cf_dmg = sum(cf_dmg, na.rm = TRUE)) %>%
        mutate(mid = obs_damage - cf_dmg)

    ## -------------------------------------------------------------------------

    diff_limits_model[[k]] <- model_sensitivity_results[[k]] %>%
         summarize(obs_damage = unique(obs_damage),
              low = obs_damage - quantile(cf_dmg, 0.975),
              high = obs_damage - quantile(cf_dmg, 0.025)) %>%
        mutate(model_code = model_sensitivity_list[k],
               mid = damage_mid$mid)
}

diff_limits_model <- bind_rows(diff_limits_model) %>%
    mutate(model_code = factor(model_code, levels = model_sensitivity_list))

## -----------------------------------------------------------------------------
## SENSITIVITY OF CUMULATIVE DAMAGES TO PRECIP TIME PERIOD
# -----------------------------------------------------------------------------
precip_sensitivity_results <- lapply(
    ppt_trend_ranges,
    function(x)
        bind_rows(lapply(list.files(path = "../processed_data/counterfactual_damages/",
                          pattern = paste0("damage_sensitivity_",
                                           x, "_region_model_"),
                          full.names = TRUE), readRDS)))

diff_limits <- vector('list', length = length(precip_sensitivity_results))

for(k in seq_along(ppt_trend_starts)){

    ## -----------------------------------------------------------------------------
    ## calculate midpoint
    quantile_breaks <- calc_bin_breaks(subset(precip_data,
                                              year >= ppt_trend_starts[k] &
                                              year <= ppt_trend_end),
                                       q_breaks,
                                       allyears = (damage_start-1):damage_end)
    quantile_trends <- calc_bin_trends(subset(precip_data,
                                              year >= ppt_trend_starts[k] &
                                              year <= ppt_trend_end),
                                       q_mid)

    detrended_precip <- detrend_precip(subset(precip_data,
                                              year >= damage_start &
                                              year <= damage_end),
                                       quantile_breaks,
                                       quantile_trends,
                                       trend_start = ppt_trend_starts[k])

    damage_mid <- calc_cf_damage(detrended_precip, damage_data,
                                 region_coeff) %>%
        summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
                  cf_dmg = sum(cf_dmg, na.rm = TRUE)) %>%
        mutate(mid = obs_damage - cf_dmg)

    ## -------------------------------------------------------------------------

    diff_limits[[k]] <- precip_sensitivity_results[[k]] %>%
         summarize(obs_damage = unique(obs_damage),
              low = obs_damage - quantile(cf_dmg, 0.975),
              high = obs_damage - quantile(cf_dmg, 0.025)) %>%
        mutate(precip_range = ppt_trend_ranges[k],
               mid = damage_mid$mid)

}
diff_limits <- bind_rows(diff_limits)

## -----------------------------------------------------------------------------
## FIGURE 3B

txt_size <- 2
pointrange_size <- 0.2
main_col <- "#279F27"
notmain_col <-  "#004D00"
main_size <- 0.5
notmain_size <- 0.3
trend_range_labels <- paste0(ppt_trend_starts,
                             c("", "", "", "\n(main)", "", "", ""))

fig_3B <- ggplot(diff_limits_model,
                                 aes(model_code)) +
    geom_linerange(aes(y = mid, ymin = low, ymax = high,
                       col = model_code, size = model_code)) +
    geom_point(size = 1.8,
                    aes(y = mid, col = model_code, size = model_code)) +
    geom_text(size = txt_size, aes(y = mid, label = paste0("  $", round(mid/1e9), "B")),
              hjust = 0) +
    geom_text(size = txt_size, aes(y = low, label = paste0("  $", round(low/1e9), "B")),
              hjust = 0) +
    geom_text(size = txt_size, aes(y = high, label = paste0("  $", round(high/1e9), "B")),
              hjust = 0) +
    scale_size_manual(values = c(main_size, rep(notmain_size, 4)), guide = "none") +
    scale_color_manual(values = c(main_col, rep(notmain_col, 4)), guide = "none") +
    scale_x_discrete(breaks = model_sensitivity_list,
                     labels = model_sensitivity_labels,
                     expand = expansion(add = c(0.3, 0.5))) +
    scale_y_continuous(name = "Billions\n(2017$)",
                       breaks = seq(30e9, 90e9,30e9),
                       labels = seq(30, 90, 30),
                       expand = c(0,0)) +
    ggtitle("Sensitivity to regression model specification:") +
    coord_cartesian(ylim = c(15e9, 110e9)) +
    theme_pub +
    theme(axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0, face = "bold"),
       plot.margin = margin(t = 4, r = 15, b = 4, l = 10, unit = "pt"))

## -----------------------------------------------------------------------------
## FIGURE 3C

fig_3C <- ggplot(diff_limits, aes(precip_range)) +
    geom_linerange(aes(y = mid, ymin = low, ymax = high,
                       col = precip_range, size = precip_range)) +
    geom_point(size = 1.8, aes(y = mid, col = precip_range)) +
    geom_text(size = txt_size, aes(y = mid, label = paste0("  $", round(mid/1e9), "B")),
              hjust = 0) +
    geom_text(size = txt_size, aes(y = low, label = paste0("  $", round(low/1e9), "B")),
              hjust = 0) +
    geom_text(size = txt_size, aes(y = high, label = paste0("  $", round(high/1e9), "B")),
              hjust = 0) +
    scale_size_manual(values = c(notmain_size, notmain_size, notmain_size, main_size,
                                 notmain_size, notmain_size, notmain_size), guide = "none") +
    scale_color_manual(values = c(notmain_col, notmain_col, notmain_col, main_col,
                                  notmain_col, notmain_col, notmain_col), guide = "none") +
    scale_x_discrete(labels = trend_range_labels,
                     expand = expansion(add = c(0.3, 0.7))) +
    scale_y_continuous(name = "Billions\n(2017$)",
                      limits = c(15e9, 110e9),
                       breaks = seq(30e9, 90e9, 30e9),
                       labels = seq(30, 90, 30),
                      expand = c(0,0)) +
    ggtitle("Sensitivity to precipitation trend starting year:") +
    theme_pub +
    theme(axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0, face = "bold"),
          plot.margin = margin(t = 4, r = 15, b = 4, l = 10, unit = "pt"))

## -----------------------------------------------------------------------------
## COMBINE FIGURE 3

pdf("./fig_3.pdf", 3.5, 3.8)
print(ggarrange(fig_3A,
                fig_3B,
                fig_3C,
                ncol = 1, nrow = 3,
                heights = c(1.35, 1, 1),
                labels = c("A", "B", "C"),
                font.label = list(size = 11, face = "plain")))
dev.off()

## -----------------------------------------------------------------------------
## SI FIG S6 -  SENSITIVITY TO DEPENDENT VARIABLE TRANSFORM
## -----------------------------------------------------------------------------

transform_sensitivity_results <- lapply(
    transform_sensitivity_list,
    function(x)
        bind_rows(lapply(list.files(path = "../processed_data/counterfactual_damages/",
                          pattern = paste0("damage_sensitivity_1928-2017_", x),
                          full.names = TRUE), readRDS)))

diff_limits_transform <- vector('list', length = length(transform_sensitivity_results))

for(k in seq_along(transform_sensitivity_list)){
    LAG <- grepl("lag", transform_sensitivity_list[k])
    QUAD <- grepl("quad", transform_sensitivity_list[k])
    POLY <- ifelse(QUAD, TRUE, FALSE)
    ## -------------------------------------------------------------------------
    ## calculate midpoint
    model_coeff <- format_coeff(readRDS(paste0("../processed_data/regression_models/",
                                       transform_sensitivity_list[k],
                                       ".Rds")),
                                n_var = 3,
                                var_codes = c("precip", "region", "season"), POLY) %>%
        dplyr::rename(coeff = Estimate) %>%
        dplyr::select(-`Std. Error`)

     ## remove season column if all NAs
    if(sum(!is.na(model_coeff$season)) == 0){
        model_coeff <- dplyr::select(model_coeff, -season)
    }
    if(sum(!is.na(model_coeff$region)) == 0){
        model_coeff<- dplyr::select(model_coeff, -region) %>%
            crossing(region = unique(damage_data$region))
    }

    damage_mid <- calc_cf_damage(base_precip_detrended, damage_data,
                                 model_coeff, LAG, QUAD) %>%
        summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
                  cf_dmg = sum(cf_dmg, na.rm = TRUE)) %>%
        mutate(mid = obs_damage - cf_dmg)

    ## -------------------------------------------------------------------------

    diff_limits_transform[[k]] <- transform_sensitivity_results[[k]] %>%
         summarize(obs_damage = unique(obs_damage),
              low = obs_damage - quantile(cf_dmg, 0.975),
              high = obs_damage - quantile(cf_dmg, 0.025)) %>%
        mutate(model_code = transform_sensitivity_list[k],
               mid = damage_mid$mid)
}

diff_limits_transform <- bind_rows(diff_limits_transform) %>%
    mutate(model_code = factor(model_code, levels = transform_sensitivity_list))

fig_S6 <- ggplot(diff_limits_transform,
                                     aes(model_code)) +
    geom_linerange(aes(y = mid, ymin = low, ymax = high,
                       col = model_code, size = model_code)) +
    geom_point(size = 1.8, aes(y = mid, col = model_code)) +
    geom_vline(xintercept = 5.65, linetype = "dashed", size = 0.3) +
    geom_text(size = txt_size, aes(y = mid, label = paste0("  $", round(mid/1e9), "B")),
              hjust = 0) +
    geom_text(size = txt_size, aes(y = low, label = paste0("  $", round(low/1e9), "B")),
              hjust = 0) +
    geom_text(size = txt_size, aes(y = high, label = paste0("  $", round(high/1e9), "B")),
              hjust = 0) +
    scale_size_manual(values = c(notmain_size, main_size, rep(notmain_size, 4)),
                      guide = "none") +
    scale_color_manual(values = c(notmain_col, main_col,
                                  rep(notmain_col, 4)),
                       guide = "none") +
    scale_x_discrete(name = "log transformations",
                     breaks = transform_sensitivity_list,
                     labels = transform_sensitivity_labels,
                     expand = expansion(add = c(0.3, 0.5))) +
    scale_y_continuous(name = "Billions (2017$)",
                       breaks = seq(30e9, 90e9, 30e9),
                       labels = seq(30, 90, 30),
                       expand = c(0,0)) +
    coord_cartesian(ylim = c(15e9, 110e9)) +
    ggtitle("Sensitivity of cumulative damages from precipitation change (Fig. 3)") +
    theme_pub +
    theme(axis.title.x = element_text(size = 7, hjust = 0.4),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(hjust = 0, face = "bold"),
           plot.margin = margin(t = 4, r = 15, b = 4, l = 10, unit = "pt"))

pdf("./fig_S6.pdf", 4, 2)
print(fig_S6)
dev.off()
