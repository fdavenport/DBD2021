library(tidyverse)
library(raster)
library(sf)
library(ggpubr)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

state_sf <- st_read(state_boundary_file) %>%
    mutate(STATENAME = NAME)

## -----------------------------------------------------------------------------
## MONTHLY PRECIPITATION TRENDS (PRISM)
quants <- c(0.5, 0.75, 0.95)

monthly_data <- bind_rows(readRDS(monthly_prism_2pt5_file))

monthly_grid_stats <- monthly_data %>%
    mutate(year = as.numeric(year)) %>%
    subset(year %in% base_period) %>%
    group_by(centroid_coord) %>%
    summarize(mean_precip = mean(monthly_precip),
              sd_precip = sd(monthly_precip))

monthly_grid_trends <- monthly_data %>%
    mutate(year = as.numeric(year)) %>%
    subset(year >= precip_start & year <= precip_end) %>%
    left_join(monthly_grid_stats) %>%
    mutate(monthly_precip_NORM = (monthly_precip - mean_precip)/sd_precip) %>%
    crossing(quant = quants) %>%
    group_by(centroid_coord, quant) %>%
    summarize(monthly_trend = rq(monthly_precip_NORM ~ year,
                                tau = unique(quant))$coefficients[2]) %>%
    ungroup() %>%
    mutate(quant = paste0("q", quant))

monthly_trend_sf <- raster_df_to_sf(
    monthly_grid_trends %>%
    spread(quant, monthly_trend)) %>%
    gather(quant, monthly_trend, contains("q0")) %>%
    mutate(trend_bin = cut(monthly_trend*10,
                           breaks = c(-Inf, seq(-0.15, 0.15, 0.05), Inf)))

monthly_trend_map <- add_map_formatting(
    ggplot(state_sf) +
    geom_sf(data = monthly_trend_sf, col = "transparent", aes(fill = trend_bin)) +
    geom_sf(color = "gray36", fill = alpha("white", 0), size = 0.2) +
    facet_grid(quant~., switch = "y", labeller=quant_labeller),
    scale_name = "trend (s.d./decade)",
    val = brewer.pal(8, "BrBG"),
    drp = FALSE,
    lab = c("-0.15", "", "-0.05", "", "0.05", "", "0.15")) +
    ggtitle("monthly precipitation\n(PRISM)") +
    theme(panel.grid.major=element_line(colour="transparent"),
          strip.text = element_text(angle = 270, colour="black", size = 7))

centroid_coords <- unique(monthly_grid_trends$centroid_coord)

## -----------------------------------------------------------------------------
## MONTHLY MAX 5-DAY TRENDS (CLIMDEX)

rx5day_data <- readRDS(rx5day_obs_file) %>%
    mutate(centroid_coord = paste0(lon, "_", lat)) %>%
    subset(centroid_coord %in% centroid_coords)

rx5day_grid_stats <- rx5day_data %>%
    mutate(year = as.numeric(year)) %>%
    subset(year %in% base_period) %>%
    group_by(centroid_coord) %>%
    summarize(mean_precip = mean(rx5day),
              sd_precip = sd(rx5day))

rx5day_trends <- rx5day_data %>%
    mutate(year = as.numeric(year)) %>%
    subset(year >= precip_start & year <= precip_end) %>%
    left_join(rx5day_grid_stats, by = c("centroid_coord")) %>%
    mutate(rx5day_NORM = (rx5day - mean_precip)/sd_precip) %>%
    crossing(quant = quants) %>%
    group_by(centroid_coord, quant) %>%
    summarize(rx5day_trend = rq(rx5day_NORM ~ year,
                                tau = unique(quant))$coefficients[2]) %>%
    ungroup() %>%
    mutate(quant = paste0("q", quant))

rx5day_trend_sf <- raster_df_to_sf(
    rx5day_trends %>%
    spread(quant, rx5day_trend)) %>%
    gather(quant, rx5day_trend, contains("q0")) %>%
    mutate(trend_bin = cut(rx5day_trend*10, breaks = c(-Inf, seq(-0.15, 0.15, 0.05), Inf)))

rx5day_map <- add_map_formatting(
    ggplot(state_sf) +
    geom_sf(data = rx5day_trend_sf, col = "transparent", aes(fill = trend_bin)) +
    geom_sf(color = "gray36", fill = alpha("white", 0), size = 0.2) +
    facet_grid(quant~., switch = "y", labeller=quant_labeller),
     scale_name = "trend (s.d./decade)",
    val = brewer.pal(8, "BrBG"),
    drp = FALSE,
    lab = c("-0.15", "", "-0.05", "", "0.05", "", "0.15")) +
    ggtitle("monthly max 5-day precipitation\n(Climdex)") +
    theme(panel.grid.major=element_line(colour="transparent"),
          strip.text = element_text(angle = 270, colour="black", size = 7))

pdf("./fig_2.pdf", 3.5, 3.5)
print(annotate_figure(ggarrange(monthly_trend_map,
                rx5day_map +
                   theme(panel.grid.major=element_line(colour="transparent"),
                         strip.text = element_text(angle = 270, colour="transparent",
                                                   size = 7)),
                ncol = 2, nrow = 1,
                labels = c("A", "B"),
                font.label = list(size = 11,
                                  face = "plain")),
                top = text_grob("Observed precipitation trends 1928-2017",
                                size = 7, face = "bold")))
dev.off()
