library(tidyverse);
library(RColorBrewer);
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
## SI FIG S2A: DISTRIBUTED LAG MODEL
lag_model <- readRDS("../processed_data/regression_models/lag_model.Rds")

lag_df <- read_boot_results("../processed_data/bootstrapped_models/lag_model_bootstrap_1.Rds",
                            n_var = 2,
                            var_codes = c("precip", "region")) %>%
    group_by(region, precip) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
                     yHi = quantile(Estimate, 0.975),
                     yMid = quantile(Estimate, 0.5)) %>%
    ungroup() %>%
    mutate(region = factor(region, levels = region_order))

lag_colors <- c("monthly_precip_NORM" = "black",
                "monthly_precip_NORM_lag1" = "dodgerblue4",
                "monthly_precip_NORM_lag2" = "dodgerblue2")
fig_S2A <- ggplot(lag_df, aes(x = region, y = yMid, col = precip)) +
    geom_hline(yintercept = 0) +
    geom_pointrange(size = 0.4, fatten = 2, aes(ymin = yLow, ymax = yHi),
                    position = position_dodge(width = 0.5)) +
    scale_y_continuous(name = expression(beta), limits = c(-.3, 2.3), expand = c(0, 0)) +
    scale_x_discrete(name = "Geographic region", labels = region_labels) +
    scale_color_manual(values = lag_colors,
                       labels = c("monthly precip. (s.d.)",
                                  "1 month lag",
                                  "2 month lag"),
                       guide = guide_legend(override.aes = list(size = .2))) +
    theme_SI +
    theme(axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          legend.title = element_blank())

pdf("./fig_S2A.pdf", 5, 2)
print(fig_S2A)
dev.off()

## -----------------------------------------------------------------------------
## READ IN REGION MODEL RESULTS

region_model_coeff <- readRDS("../processed_data/regression_models/region_model.Rds") %>%
    format_coeff(., n_var = 2, var_codes = c("precip", "region")) %>%
    dplyr::select(region, yMid = Estimate)
region_df <- read_boot_results("../processed_data/bootstrapped_models/region_model_bootstrap_1.Rds",
                                       n_var = 2,
                                       var_codes = c("precip", "region")) %>%
    group_by(region) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
              yHi = quantile(Estimate, 0.975)) %>%
    left_join(region_model_coeff) %>%
    mutate(region = factor(region, levels = region_order))

region_ribbon_df <- state_order %>%
    left_join(region_df) %>%
    mutate(abb = factor(abb, levels = state_order$abb),
           region = factor(region, levels = region_order)) %>%
    na.omit()
## -----------------------------------------------------------------------------
## SI FIG S2B: TEMPORAL BLOCK BOOTSTRAP

region_tbb_df <- read_boot_results("../processed_data/bootstrapped_models/region_tbb_bootstrap_1.Rds",
                                   n_var = 2,
                                   var_codes = c("precip", "region")) %>%
    group_by(region) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
              yHi = quantile(Estimate, 0.975)) %>%
    left_join(region_model_coeff) %>%
    mutate(region = factor(region,
                           levels = region_order))

fig_S2B <- ggplot(region_df, aes(x = region, y = yMid)) +
    geom_pointrange(size = 0.4, fatten = 2, aes(ymin = yLow, ymax = yHi,
                                                col = "black")) +
    geom_pointrange(size = 0.4, fatten = 2, data = region_tbb_df,
                    aes(x = as.numeric(region) + 0.2, ymin = yLow, ymax = yHi,
                    col = "green4")) +
    scale_y_continuous(name = expression(beta), limits = c(0, 2.3), expand = c(0, 0)) +
    scale_x_discrete(name = "Geographic region", labels = region_labels) +
    scale_color_manual(name = "",
                       values = c("black", "green4"),
                       labels = c("bootstrapping states",
                                  "temporal block bootstrap (3-year block)"),
                       guide = guide_legend(override.aes = list(size = .2))) +
    theme_SI +
    theme(axis.ticks.x = element_blank(),
          legend.position = c(0.7, .2),
          legend.background = element_blank())

pdf("./fig_S2B.pdf", 4, 2)
print(fig_S2B)
dev.off()

## -----------------------------------------------------------------------------
## SI FIG S2C: DECADE MODEL

decade_model <- readRDS("../processed_data/regression_models/decade_model.Rds")

decade_df <- read_boot_results("../processed_data/bootstrapped_models/decade_model_bootstrap_1.Rds",
                                       n_var = 2,
                                       var_codes = c("precip", "decade")) %>%
    group_by(decade) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
              yHi = quantile(Estimate, 0.975)) %>%
    merge(format_coeff(decade_model, n_var = 2, var_codes = c("precip", "decade")) %>%
          dplyr::select(decade, yMid = Estimate))

fig_S2C <- ggplot(decade_df, aes(x = decade, y = yMid)) +
    geom_pointrange(size = 0.3, fatten = 1.5, aes(ymin = yLow, ymax = yHi)) +
    scale_y_continuous(name = expression(beta), limits = c(0.5, 1.8), expand = c(0, 0)) +
    scale_x_discrete(name = "Time Period") +
    theme_pub +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(size = 7, angle = 0),
          plot.margin = margin(t = 10, b = 2, l = 5, r = 3, unit = "pt"))

pdf("./fig_S2C.pdf", 2.5, 2)
print(fig_S2C)
dev.off()

## -----------------------------------------------------------------------------
## SI FIG S2D: STATE MODEL
state_coeff <- readRDS("../processed_data/regression_models/state_model.Rds") %>%
    format_coeff(., n_var = 2, var_codes = c("precip", "STATENAME")) %>%
    dplyr::rename(yMid = Estimate)

state_df <- read_boot_results("../processed_data/bootstrapped_models/state_model_bootstrap_1.Rds",
                              n_var = 2,
                              var_codes = c("precip", "STATENAME")) %>%
    group_by(STATENAME, precip) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
                     yHi = quantile(Estimate, 0.975)) %>%
    ungroup() %>%
    mutate(state_abb = state.abb[match(STATENAME, state.name)],
           state_abb = factor(state_abb, levels = state_order$abb),
           region = factor(findNCEIregion(STATENAME),
                           levels = region_order)) %>%
    left_join(state_coeff)

fig_S2D <- ggplot(subset(state_df, precip == "monthly_precip_NORM")) +
    geom_ribbon(data = region_ribbon_df, alpha = 0.2,
                aes(ymin = yLow, ymax = yHi, y = NULL,
                    x = as.numeric(abb), fill = region)) +
    geom_line(data = region_ribbon_df, col = "black",
              aes(x = as.numeric(abb), y = yMid,
                  group = region), size = 0.5,
              linetype = "dotted") +
    geom_pointrange(size = 0.3, fatten = 1, col = "black",
                    aes(x = as.numeric(state_abb), y = yMid,
                        ymin = yLow, ymax = yHi)) +
    geom_text(size = 2.5, aes(x = as.numeric(state_abb), y = yHi + .15,
                            label = state_abb), angle = 90) +
    scale_y_continuous(name = expression(beta), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0,0),
                       breaks = region_spacing,
                       labels = region_label_long) +
    scale_fill_manual(values = rep("gray36", 10))+
    coord_cartesian(ylim = c(0, 2.9), xlim = c(0, 49)) +
    theme_SI +
    theme(axis.text.x = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 8))


pdf("./fig_S2D.pdf", 6.5, 2.5)
print(fig_S2D)
dev.off()

