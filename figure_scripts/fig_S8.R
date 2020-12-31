library(tidyverse)
library(ggpubr)
library(lfe)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

damage_data <- readRDS(state_panel_file)

## -----------------------------------------------------------------------------
## calculate detrended precip

precip_data <- readRDS(precip_std_file)

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

precip_demeaned <- demeanlist(dplyr::select(damage_data, monthly_precip_NORM),
                          list(factor(damage_data$STATEYEAR), factor(damage_data$STATEMONTH)))

fill_col <- c("observed" = "gray60",
              "difference" = "dodgerblue4")
fig_S8A <- ggplot(precip_demeaned, aes(value)) +
    geom_histogram(data = precip_demeaned, alpha = 0.5, binwidth = 0.2,
                   aes(x = monthly_precip_NORM, fill = "observed")) +
    geom_histogram(data = subset(detrended_precip, year >= damage_start),
                   alpha = 0.5, binwidth = 0.2,
                   aes(x = precip_diff, fill = "difference")) +
    scale_fill_manual(values = fill_col,
                      breaks = c("observed", "difference"),
                      labels = c("observed precipitation\n(after FE removed)",
                                 "change in precipitation\n(counterfactual)")) +
    scale_y_continuous(name = "Number of state-months") +
    scale_x_continuous(name = "monthly precipitation anomaly (std. dev.)",
                       limits = c(-2.5, 5.5)) +
    theme_SI +
    theme(legend.title = element_blank(),
          legend.position = c(0.8, 0.8),
          legend.key = element_rect(size = 3, color = "white"),
          legend.key.size = unit(1, 'lines')) +
    guides(fill = guide_legend(override.aes = list(size = 2)))


fig_b_data <- detrended_precip %>%
    dplyr::select(STATENAME, monthly_precip_NORM, monthly_precip_NORM_detrend) %>%
    gather(key, value, c("monthly_precip_NORM", "monthly_precip_NORM_detrend"))

dist_col <- c("monthly_precip_NORM" = "gray60",
              "monthly_precip_NORM_detrend" = "#2171B5")
outline_col <- c("monthly_precip_NORM" = "black",
                 "monthly_precip_NORM_detrend" = "#2171B5")

fig_S8B <- ggplot(fig_b_data, aes(value)) +
    geom_density(alpha = 0.5, col = "transparent", aes(fill = key)) +
    stat_density(aes(col = key), geom="line", position = "identity") +
    scale_fill_manual(values = dist_col,
                      guide = "none") +
    scale_colour_manual(values = outline_col,
                        labels = c("observed precipitation", "detrended precipitation")) +
    scale_x_continuous(name = "monthly precipitation anomaly (std. dev.)",
                       limits = c(-3, 5.5)) +
    theme_SI+
    theme(legend.title = element_blank(),
          legend.position = c(0.8, 0.8),
           legend.key.size = unit(1, 'lines'))


pdf("./fig_S8.pdf", 3.5, 4.5)
print(ggarrange(fig_S8A, fig_S8B,
                nrow = 2, ncol = 1,
                heights = c(2.25, 2),
                labels = c("A", "B"),
                align = "v",
                font.label = list(size = 11, face = "plain")))
dev.off()

