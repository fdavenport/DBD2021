library(tidyverse)
library(raster)
library(sf)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
data <- readRDS(state_panel_file)

## -----------------------------------------------------------------------------
## FIGURE 1A: MAP OF STATE  DAMAGE TRENDS
statebounds <- shapefile(state_boundary_file)
statebounds$STATENAME <- statebounds$NAME
statebounds$region <- findNCEIregion(statebounds$STATENAME)
state_sf <- st_as_sf(statebounds)

state_damage_trends <- data %>%
    group_by(STATENAME) %>%
    summarize(damage_trend = lm(damage_value ~ year)$coefficients[2])

trend_sf <- state_sf %>%
    merge(state_damage_trends) %>%
    mutate(d_trend_bin = cut(damage_trend,
                             breaks = c(-Inf, 0, .01, .02, .04, .08, Inf)))

## -----------------------------------------------------------------------------
## region spatial data

region_sf <- st_as_sf(aggregate(statebounds, by = "region")) %>%
    subset(region != "otherNCEI")

region_label_df <- data.frame(region = region_sf$region,
                              st_coordinates(st_centroid(region_sf))) %>%
    mutate(xadj = c(8.5, 9, -3, -12, 4, 4.5, -3, 2.5, -8),
           X_new = X + xadj,
           yadj = c(8.5, -3, 6, 0, -7, -4, -8, 6, -3),
           Y_new = Y + yadj,
           region_id = c("C", "NE", "NR",
                         "NW", "S", "SE",
                         "SW", "UM","W"))

region_label_sf <- st_as_sf(region_label_df,
                            coords = c("X_new", "Y_new"),
                            crs = crs(region_sf))

line_sf <- st_cast(st_linestring(x = matrix(c(-87.16669, 38.71032,
                                      -78.66669, 47.21032),
                                      2, 2), "XY"),
                   to= "MULTILINESTRING",
                    crs = st_crs(region_sf))

line_sf <- data.frame(id = c("Central", "Central"),
                      X = c(-87.16669, -80.5),
                      Y = c(38.71032, 45.5)) %>%
    st_as_sf(coords = c("X", "Y")) %>%
    group_by(id) %>%
    summarize() %>%
    st_cast("LINESTRING") %>%
    st_set_crs(crs(region_label_sf))

## -----------------------------------------------------------------------------

fig_1A <- ggplot(trend_sf, aes(fill = d_trend_bin)) +
    geom_sf(color = "gray45", size = 0.1) +
    geom_sf(data = region_sf, col = "gray20", fill = NA, size = 0.3) +
    geom_sf_text(data = region_label_sf, size = 2,
                 aes(label = region_id, fill = NULL)) +
    geom_sf(data = line_sf, size = 0.3, lineend = "round",
            aes(fill = NULL)) +
    scale_fill_manual(name =  "\n\nln(norm. damages)\n per year",
                      values = brewer.pal(6, "Oranges"),
                      labels = c("0", "0.01", "0.02", "0.04", "0.08"),
                      drop = FALSE) +
    coord_sf(crs = st_crs(2163), xlim = c(-2650000, 2800000),
             ylim = c(-2450000,805000)) +
    ggtitle(paste0("Trends in monthly flood damages\n(", damage_start, "-2017)")) +
    map_theme_pub +
    theme(legend.position = "right",
          panel.grid.major=element_line(colour="transparent"),
          legend.title = element_text(angle = 90, vjust = 0.5,
                                      hjust = 0.5, size = 6.5),
          legend.text = element_text(size = 6.5),
          plot.title = element_text(size = 7, hjust = 0.5),
          plot.margin = margin(t = -4, r = -2, b = -6, l = -2, unit = "pt")) +
    guides(fill = guide_coloursteps(barwidth = 0.2,
                                    barheight = 4,
                                    ticks = TRUE,
                                    ticks.colour = "gray36",
                                    title.position = "left",
                                    frame.colour = "gray36",
                                    frame.linewidth = .1))

## -----------------------------------------------------------------------------
## REGRESSION MODEL RESULTS
## -----------------------------------------------------------------------------
x_mean <- mean(data$monthly_precip_NORM)
y_mean <- mean(data$damage_value)

## READ IN POOLED LINEAR MODEL
linear_model <- readRDS("../processed_data/regression_models/linear.Rds")
lm_slope <- linear_model$coefficients[1]
lm_int <- y_mean-(lm_slope*x_mean)
linear_boot_coeff <- read_boot_results("../processed_data/bootstrapped_models/linear_bootstrap_1.Rds",
                                       n_var = 1,
                                       var_codes = "precip")
lm_slope_low <- quantile(linear_boot_coeff$Estimate, 0.025)
lm_int_low <- y_mean-(lm_slope_low*x_mean)
lm_slope_hi <- quantile(linear_boot_coeff$Estimate, 0.975)
lm_int_hi <- y_mean-(lm_slope_hi*x_mean)

x <- seq(-2.5, 5, 0.5)
linear_df <- data.frame(x) %>%
    mutate(yMid = x*lm_slope + lm_int,
           yLow = x*quantile(linear_boot_coeff$Estimate, 0.025) + lm_int_low,
           yHi = x*quantile(linear_boot_coeff$Estimate, 0.975) + lm_int_hi)

## READ IN BINNED MODEL
model_nl <- readRDS("../processed_data/regression_models/binned_model.Rds")
nl_boot_results <- readRDS("../processed_data/bootstrapped_models/binned_bootstrap_1.Rds")
nl_boot_coeff <- do.call(bind_rows, lapply(nl_boot_results, function(x) {
    coeff_offset <- y_mean - sum(x[4:5])/2
    as.data.frame(x) %>% rownames_to_column(var = "precipBin") %>%
            mutate(ref_estimate = Estimate + coeff_offset)}))

nl_df <-  nl_boot_coeff %>%
    mutate(precipBin = gsub("precip_binstd\\(", "", precipBin)) %>%
    group_by(precipBin) %>%
    summarize(yMid = quantile(ref_estimate, 0.5),
           yLow = quantile(ref_estimate, 0.025),
           yHi = quantile(ref_estimate, 0.975)) %>%
    ungroup() %>%
    separate(precipBin, sep = ",", into = c("xStart", "xEnd")) %>%
    mutate(xStart = as.numeric(xStart),
           xEnd = as.numeric(gsub("]", "", xEnd)),
           xMid = (xEnd + xStart)/2,
           xRibbon = ifelse(xMid == -1.75, -2, xMid)) %>%
    arrange(xStart)

## READ IN QUADRATIC MODEL
quad_coeff <- format_coeff(readRDS("../processed_data/regression_models/quad_model.Rds"),
                           n_var = 1, var_codes = "precip", POLY = TRUE) %>%
    dplyr::select(-`Std. Error`) %>%
    pivot_wider(names_from = poly_order, names_prefix = "poly", values_from = Estimate)
quad_coeff_boot <- read_boot_results(
    "../processed_data/bootstrapped_models/quad_model_bootstrap_1.Rds",
    n_var = 1, var_codes = "precip", POLY = TRUE) %>%
    dplyr::select(-`Std. Error`) %>%
    pivot_wider(names_from = poly_order, names_prefix = "poly", values_from = Estimate)
quad_df <- crossing(x, quad_coeff_boot) %>%
    mutate(y = poly1*x + poly2*x^2) %>%
    group_by(x) %>%
    summarize(yLow = quantile(y, 0.025) + lm_int,
              yHi = quantile(y, 0.975) + lm_int) %>%
    ungroup() %>%
    mutate(poly1 = quad_coeff$poly1,
           poly2 = quad_coeff$poly2,
           x = as.numeric(as.character(x)),
           yMid = poly1*x + poly2*x^2 + lm_int)

## -----------------------------------------------------------------------------
## FIGURE 1B: PANEL MODELS
## -----------------------------------------------------------------------------
ribbon_col <- c("linear" = "#6BAED6",
                "binned" = "brown3",
                "quadratic" = "gray40")
line_col <- c("linear" = "#08519C",
              "binned" = "brown3",
              "quadratic" = "gray30")
fig_1B <- ggplot(linear_df, aes(x = x, y = yMid)) +
    geom_ribbon(data = quad_df, alpha = 0.4,
                aes(x = x, ymin = yLow, ymax = yHi, fill = "quadratic")) +
    geom_line(data = quad_df, size = 0.4, aes(x, yMid, col = "quadratic")) +
    geom_ribbon(data = nl_df, alpha = 0.4,
                aes(x = xRibbon, ymin = yLow, ymax = yHi, fill = "binned")) +
    geom_line(data = nl_df, size = 0.4, aes(xMid, yMid, col = "binned")) +
    geom_point(data = nl_df, size = 0.6, col = "brown3", aes(xMid, yMid)) +
    geom_ribbon(alpha = 0.6, aes(ymin = yLow, ymax = yHi, fill = "linear")) +
    geom_line(size = 0.4, aes(col = "linear")) +
    scale_fill_manual(values = ribbon_col, guide = "none") +
    scale_color_manual(values = line_col,
                       breaks = c("linear", "quadratic", "binned")) +
    scale_y_continuous(name = "ln(normalized damages)", expand = c(0, 0),
                       breaks = seq(-4, 8, 2), labels = seq(-4, 8, 2)) +
    scale_x_continuous(name = "Monthly precipitation anomaly (s.d.)",
                       breaks = seq(-2.5, 5, 2.5),
                       labels = c("-2.5", "0", "2.5", "5"),
                       expand = c(0, 0)) +
    coord_cartesian(xlim = c(-2, 5)) +
    ggtitle("Flood damages vs. precipitation\n(state-month panel regression)") +
    theme_pub +
    theme(legend.position = c(0.3, 0.85),
          legend.title = element_blank(),
          legend.key.size = unit(7, "pt"),
          axis.line.x = element_line(size = 0.2),
          plot.margin = margin(l = 0.5, r = 5, unit = "pt"),
          axis.title.x = element_blank())

## -----------------------------------------------------------------------------
## HISTOGRAMS
## -----------------------------------------------------------------------------
precip_hist <- ggplot(data, aes(x = monthly_precip_NORM)) +
    geom_histogram(alpha = 0.6, binwidth = 0.25, fill = "#6BAED6") +
scale_x_continuous(name = "Precipitation anomaly (s.d.)",
                   breaks = seq(-2.5, 5, 2.5),
                   labels = c("-2.5", "0", "2.5", "5"),
                   expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("dist. of monthly precipitation:") +
    theme_pub +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(size = 0.2),
          axis.ticks.y = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(b = 0, unit = "pt"),
          plot.title = element_text(size = 6.5, hjust = -0.7, vjust = 0,
                                    margin=margin(b = 2, unit = "pt"))) +
    coord_cartesian(ylim = c(0, 2000), xlim = c(-2, 5))

damage_hist <- ggplot(subset(data, damage_binary == 1), aes(x = monthly_precip_NORM)) +
    geom_histogram(alpha = 0.6, binwidth = 0.25, fill = "gray51") +
    scale_x_continuous(name = "Precipitation anomaly (s.d.)",
                          breaks = seq(-2.5, 5, 2.5),
                       labels = c("-2.5", "0", "2.5", "5"),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("dist. of months with damage:") +
    theme_pub +
    theme(axis.title = element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_line(size = 0.2),
          axis.ticks.y = element_blank(),
          axis.text = element_blank(),
          plot.margin = margin(b = 0, unit = "pt"),
          plot.title = element_text(size = 6.5, hjust = -0.7, vjust = 0,
                                    margin=margin(b = 2, unit = "pt"))) +
    coord_cartesian(ylim = c(0, 670), xlim = c(-2, 5))

dollar_hist <- ggplot(data, aes(x = monthly_precip_NORM, weights = totaldmg_adj2017)) +
    geom_histogram(alpha = 0.6, binwidth = 0.3, fill = "black") +
    scale_x_continuous(name = "Monthly precipitation anomaly (s.d.)",
                          breaks = seq(-2.5, 5, 2.5),
                       labels = c("-2.5", "0", "2.5", "5"),
                       expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggtitle("dist. of damages ($):") +
    theme_pub +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 7, margin = margin(t = 0, unit = "pt")),
          axis.line.y = element_blank(),
          axis.line.x = element_line(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(b = 0, unit = "pt"),
          plot.title = element_text(size = 6.5, hjust = -0.2, vjust = 0,
                                    margin=margin(b = -3, unit = "pt"))) +
    coord_cartesian(xlim = c(-2, 5), ylim = c(0,50e9))


## -----------------------------------------------------------------------------
## FIGURE 1C: REGIONAL MODELS
## -----------------------------------------------------------------------------
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

region_5day_coeff <- readRDS("../processed_data/regression_models/region_5day_model.Rds") %>%
    format_coeff(., n_var = 2, var_codes = c("precip", "region")) %>%
    dplyr::select(region, yMid = Estimate)

region_5day_df <- read_boot_results("../processed_data/bootstrapped_models/region_5day_model_bootstrap_1.Rds",
                                   n_var = 2,
                                   var_codes = c("precip", "region")) %>%
    group_by(region) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
              yHi = quantile(Estimate, 0.975)) %>%
    left_join(region_5day_coeff) %>%
    mutate(region = factor(region,
                           levels = region_order))

## -----------------------------------------------------------------------------
## check statistical significance of regional models
region_diff <- read_boot_results("../processed_data/bootstrapped_models/region_model_diff.Rds",
                               n_var = 2,
                       var_codes = c("precip", "region")) %>%
    group_by(region) %>%
    summarize(ecdf_0 = ecdf(damage_value)(0)) %>%
    mutate(p_diff = ifelse(ecdf_0 < 0.5, ecdf_0*2, (1-ecdf_0)*2),
           symbol = ifelse(p_diff < 0.05, 19, 21)) %>%
    mutate(region = factor(region,
                           levels = region_order))

print(region_diff)

region_5day_diff <- read_boot_results("../processed_data/bootstrapped_models/region_5day_model_diff.Rds",
                                      n_var = 2,
                                      var_codes = c("precip", "region")) %>%
    group_by(region) %>%
    summarize(ecdf_0 = ecdf(damage_value)(0)) %>%
    mutate(p_diff = ifelse(ecdf_0 < 0.5, ecdf_0*2, (1-ecdf_0)*2),
           symbol = ifelse(p_diff < 0.05, 19, 21)) %>%
    mutate(region = factor(region,
                           levels = region_order))

sig_legend <- get_legend(ggplot(region_diff, aes(region, p_diff)) +
    geom_point(size = 0.6, aes(color = as.factor(symbol))) +
    scale_color_manual(guide = "legend",
                       values = c("gray50", "black"),
                       labels = c("", " p < 0.05")) +
    theme_pub +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 6.5),
          legend.key.size = unit(1, "pt"),
          legend.spacing.x = unit(1.5, "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt")))

print(region_5day_diff)

region_df <- region_df %>% left_join(region_diff[c("region", "symbol")], by = "region")
region_5day_df <- region_5day_df %>% left_join(region_5day_diff[c("region", "symbol")], by = "region")

fig_1C <- ggplot(region_df, aes(x = region, y = yMid)) +
    geom_hline(yintercept = lm_slope, size = 0.3, linetype = "dashed") +
    geom_linerange(size = 0.2, aes(ymin = yLow, ymax = yHi, col = "black")) +
    geom_point(size = 0.6, fill = "white", aes(col = "black", shape = symbol)) +
    geom_linerange(size = 0.2, data = region_5day_df,
                    aes(x = as.numeric(region) + 0.25, ymin = yLow, ymax = yHi,
                        col = "gray50")) +
    geom_point(size = 0.6, fill = "white", data = region_5day_df,
                   aes(x = as.numeric(region) + 0.25,
                       col = "gray50", shape = symbol)) +
    scale_shape_identity(guide = "none") +
    scale_y_continuous(name = expression(beta), limits = c(0, 2.5), expand = c(0, 0),
                       breaks = c(0, 1, 2)) +
    scale_x_discrete(name = "Geographic region", labels = region_labels) +
    scale_color_manual(name = "",
                       values = c("black", "gray50"),
                       labels = c("monthly  ",
                                  "monthly max 5-day"),
                       guide = guide_legend(override.aes = list(size = .35))) +
    ggtitle(expression("Regional variation in"~beta)) +
    theme_pub +
    theme(axis.ticks.x = element_blank(),
          legend.position = c(0.5, .95),
          axis.title.x = element_blank(),
          axis.line.y = element_line(size = 0.3),
          legend.title = element_blank(),
          legend.text = element_text(size = 6.5),
          legend.direction = "horizontal",
          legend.key.size = unit(1, "pt"),
          legend.background = element_blank(),
          plot.title = element_text(size = 7, hjust = 0.5),
          plot.margin = margin(r = 0, l = 3, b = 6, unit = "pt"),
          axis.text.x = element_text(angle = 90),
          axis.title.y = element_text(angle = 0,
                                      vjust = 0.5, margin = margin(r = 0, unit = "pt")))

## -----------------------------------------------------------------------------
## FIGURE 1D: SEASONAL MODELS
## -----------------------------------------------------------------------------
region_ribbon_df <- state_order %>%
    left_join(region_df) %>%
    mutate(abb = factor(abb, levels = state_order$abb),
           region = factor(region, levels = region_order)) %>%
    na.omit()

season_coeff <- readRDS("../processed_data/regression_models/season_model.Rds") %>%
    format_coeff(., n_var = 3, var_codes = c("precip", "region", "season")) %>%
    dplyr::rename(yMid = Estimate)

season_df <- read_boot_results("../processed_data/bootstrapped_models/season_model_bootstrap_1.Rds",
                               n_var = 3,
                               var_codes = c("precip", "region", "season")) %>%
    group_by(region, precip, season) %>%
    dplyr::summarize(yLow = quantile(Estimate, 0.025),
                     yHi = quantile(Estimate, 0.975)) %>%
    ungroup() %>%
    mutate(region = factor(region, levels = region_order),
           season = factor(season, levels = season_order)) %>%
    merge(season_coeff)


fig_1D <- ggplot(na.omit(season_df)) +
    geom_rect(data = region_ribbon_df, alpha = 0.05, fill = "gray36",
              aes(ymin = yLow, ymax = yHi, y = NULL,
                  xmin = as.numeric(region) - 0.4,
                  xmax = as.numeric(region) + 0.4)) +
    geom_segment(data = region_ribbon_df, col = "black", linetype = "dotted",
                 aes(x = as.numeric(region) - 0.4,
                     xend = as.numeric(region) + 0.4,
                     y = yMid, yend = yMid),
                 size = 0.3) +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_pointrange(size = 0.2, fatten = 0.4,
                    aes(x = as.numeric(region),
                        y = yMid, ymin = yLow, ymax = yHi,
                        col = season),
                    position = position_dodge2(width = 0.7)) +
    scale_color_manual(values = season_col,
                       guide = guide_legend(override.aes = list(size = 0.05))) +
    scale_x_continuous(breaks = 1:9, labels = region_labels,
                       expand = expansion(add = c(.2, .2))) +
    coord_cartesian(ylim = c(-0.5, 3)) +
    scale_y_continuous(name = expression(beta), expand = c(0, 0)) +
    theme_pub +
    ggtitle(expression("Seasonal variation in"~beta)) +
    theme(axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(size = 0.3),
          legend.title = element_blank(),
          legend.text = element_text(size = 6.5),
          legend.position = c(0.6, 0.95),
          legend.direction = "horizontal",
          legend.key.size = unit(1, "pt"),
          legend.background = element_blank(),
          plot.title = element_text(size = 7, hjust = 0.5, face = "bold"),
          plot.margin = margin(r = 0, l = 3, b = 6, unit = "pt"),
          axis.text.x = element_text(angle = 90,
                                     margin = margin(t = -3,
                                                     unit = "pt"),
                                     hjust = 1,
                                     vjust = 0.5),
          axis.title.y = element_text(angle = 0,
                                      vjust = 0.5, margin = margin(r = 0, unit = "pt")))


## -----------------------------------------------------------------------------
## COMBINED FIGURE 1

pdf("./fig_1.pdf", 3.5, 4.15)
print(ggarrange(
    ggarrange(NULL,
              fig_1A,
              NULL, ncol = 3, nrow = 1, widths = c(0.35, 1, 0.15),
              labels = c("", "A", ""), font.label = list(size = 11, face = "plain"),
              hjust = 1),
    NULL,
    ggarrange(
        ggarrange(fig_1B,
                  precip_hist, damage_hist, dollar_hist,
                  align = "v",
                  ncol = 1, nrow = 4, heights = c(.8, 0.15, .15, .26),
                  labels = c("B", "", "", ""), font.label = list(size = 11, face = "plain"),
                  vjust = 1.1),
        ggarrange( ggdraw(fig_1C) +
                   draw_plot(sig_legend, .25, .25, .1, .08),
                  NULL,
                  fig_1D,
                  ncol = 1, nrow = 3, heights = c(1, .05, 1),
                  labels = c("C","",  "D"),
                  font.label = list(size = 11, face = "plain"),
                  vjust = c(1.1, 0, 0.5, 0, 0.5), hjust = -1),
        ncol = 2, nrow = 1, widths = c(0.98, 1)),
    ncol = 1, nrow = 3, heights = c(1.25, .1,  2.8)))
dev.off()

