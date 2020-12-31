library(tidyverse)
library(raster)
library(sf)
library(lfe)
library(ggpubr)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
## INPUT DATA
data <- readRDS(state_panel_file)

## WATERSHED PRECIPITATION TIME SERIES
watershed_dat <- readRDS(watershed_precip_file)
watershed_stats <- watershed_dat %>%
    subset(year %in% base_period) %>%
    group_by(NAME) %>%
    summarize(mean_precip = mean(monthly_precip),
              sd_precip = sd(monthly_precip))
watershed_dat <- watershed_dat %>%
    left_join(watershed_stats) %>%
    mutate(watershed_precip_NORM = (monthly_precip - mean_precip)/sd_precip) %>%
    dplyr::select(-monthly_precip, -mean_precip, -sd_precip) %>%
    spread(NAME, watershed_precip_NORM)

MO_dat <- subset(data, STATENAME == "Missouri") %>%
    left_join(watershed_dat)

## -----------------------------------------------------------------------------
## CREATE MAP
statebounds <- shapefile(state_boundary_file)
state_sf <- st_as_sf(statebounds)
MO_bounds <- subset(statebounds, NAME == "Missouri")

## Mississippi river boundary
MIR_bounds <- shapefile(um_basin_file)
## Missouri river boundary
MOR_bounds <- shapefile(mo_basin_file)

MO_sf <- st_as_sf(MO_bounds)
MOR_sf <- st_as_sf(MOR_bounds)
MIR_sf <- st_as_sf(MIR_bounds)
MIR_noMO_sf <- st_as_sf(MIR_bounds - MO_bounds)
MOR_noMO_sf <- st_as_sf(MOR_bounds - MO_bounds)

basin_col <- c("#a6cee3", "#1f78b4", "#807DBA")
names(basin_col) <- c("MOR", "MO", "MIR")

## -----------------------------------------------------------------------------
## SI FIG S3A:

fig_S3A <- ggplot(state_sf) +
    geom_sf(fill = "transparent", col = "gray36") +
    geom_sf(data = MO_sf, aes(fill = "MO")) +
    geom_sf(data = MOR_noMO_sf, aes(fill = "MOR"), alpha = 0.8) +
    geom_sf(data = MIR_noMO_sf, aes(fill = "MIR"), alpha = 0.8) +
    geom_sf(data = MIR_sf, fill = "transparent", col = "gray28") +
    geom_sf(data = MOR_sf, fill = "transparent", col = "gray28") +
    geom_text(size = 2.5, label = "Upper Mississippi\nBasin",
              x = 1500000, y = 500000) +
    geom_text(size = 2.5, label = "Missouri Basin",
              x = -500000, y = 700000) +
    geom_segment(size = 0.3, x = -20000, xend = -145000,
                 y = 5000, yend = 600000) +
    geom_segment(size = 0.3, x = 650000, xend = 1055000,
                 y = 5000, yend = 495000) +
    scale_fill_manual(values = basin_col, guide = "none") +
    coord_sf(crs = st_crs(2163), xlim = c(-2300000, 2500000),
             ylim = c(-2100000,950000)) +
    map_theme_pub +
    theme(panel.grid.major=element_line(colour="transparent"))

## -----------------------------------------------------------------------------
## FIT MODEL FOR MISSOURI

MO_model <- felm(damage_value ~ monthly_precip_NORM + Mississippi_noMO +
                Missouri_noMO | STATEMONTH + STATEYEAR, data = MO_dat)

MO_coeff <- as.data.frame(summary(MO_model)$coefficients) %>%
    rownames_to_column(var = "area_name") %>%
    dplyr::select(area_name, coeff = Estimate) %>%
    mutate(area_name = ifelse(area_name == "monthly_precip_NORM", "Missouri", area_name))

## -----------------------------------------------------------------------------
## BOOTSTRAP MODEL
years <- as.character(unique(MO_dat$year))
data_split <- split(MO_dat, f = MO_dat[,"year"])
sample_list <- mbb(years, block_size = 3, nboot = 1000)
boot_results <- vector('list', length = 1000)
for(i in 1:1000){
    if(i%%100 == 0) print(i)
    current_data <- do.call(bind_rows, data_split[sample_list[[i]]])
    current_model <- felm(damage_value ~ monthly_precip_NORM + Mississippi_noMO +
                Missouri_noMO | STATEYEAR + STATEMONTH, data = current_data)
    boot_results[[i]] <- summary(current_model)$coefficients
}
names(boot_results) <- paste("boot", 1:1000)

pvals <- bind_rows(lapply(boot_results, function(x) as.data.frame(x) %>%
                                 rownames_to_column(var = "vars"))) %>%
    dplyr::select(vars, Estimate) %>%
    group_by(vars) %>%
    summarize(p_ecdf = ecdf(Estimate)(0)*2) ## two-sided
print(pvals)

coeff_ranges <- bind_rows(lapply(boot_results, function(x) as.data.frame(x) %>%
                                 rownames_to_column(var = "vars"))) %>%
    dplyr::select(vars, Estimate) %>%
    group_by(vars) %>%
    summarize(mid = median(Estimate),
              high = quantile(Estimate, 0.975),
              low = quantile(Estimate, 0.025)) %>%
    mutate(vars = ifelse(vars == "monthly_precip_NORM", "Missouri", vars),
           vars = factor(vars, levels = c("Missouri", "Missouri_noMO",
                                    "Mississippi_noMO")))

## -----------------------------------------------------------------------------
## SI FIG S3B:

fig_S3B <- ggplot(coeff_ranges, aes(x = vars, y = mid)) +
    geom_pointrange(aes(ymin = low, ymax = high, col = vars)) +
    geom_hline(yintercept = 0) +
    ylab(expression(beta)) +
    ggtitle("Estimated effect of precipitation on Missouri\n flood damages (regression coefficients)") +
    scale_x_discrete(labels = c("Missouri\nprecipitation",
                                "Missouri Basin\nprecipitation",
                                "Upper Mississippi\n precipitation")) +
    scale_color_manual(values = c("#1f78b4", "#a6cee3",  "#807DBA"),
                       guide = "none") +

    theme_SI +
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())

## -----------------------------------------------------------------------------
## COUNTERFACTUAL DAMAGES

precip_data <- readRDS(precip_std_file) %>%
    subset(STATENAME == "Missouri") %>%
    ungroup() %>%
    dplyr::select(year, month, Missouri = monthly_precip_NORM) %>%
    left_join(watershed_dat) %>%
    gather(STATENAME, monthly_precip_NORM, 3:7) %>%
    subset(year >= precip_start & year <= precip_end)

quantile_breaks <- calc_bin_breaks(precip_data, q_breaks,
                                   allyears = (damage_start-1):damage_end)
quantile_trends <- calc_bin_trends(precip_data, q_mid)
detrended_precip <- detrend_precip(precip_data,
                                   quantile_breaks,
                                   quantile_trends,
                                   trend_start = precip_start)
## region model CF damages
region_coeff <- format_coeff(readRDS("../processed_data/regression_models/region_model.Rds"),
                              n_var = 2,
                             var_codes = c("precip", "region")) %>%
    dplyr::select(region, precip, coeff = Estimate)

region_cf_dmg <- subset(detrended_precip, STATENAME == "Missouri") %>%
    dplyr::select(STATENAME, year, month, precip_diff) %>%
    inner_join(data) %>%
    left_join(region_coeff) %>%
    mutate(beta_prod = coeff*(-precip_diff),
           change_factor = exp(beta_prod),
           cf_dmg = totaldmg_adj2017*change_factor) %>%
    summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
              cf_dmg = sum(cf_dmg, na.rm = TRUE)) %>%
    mutate(diff = obs_damage - cf_dmg,
           model = "region")

## MO basin model CF damages
basin_cf_dmg <- detrended_precip %>%
    mutate(area_name = STATENAME,
           STATENAME = "Missouri") %>%
    dplyr::select(STATENAME, area_name, year, month, precip_diff) %>%
    inner_join(data) %>%
    inner_join(MO_coeff) %>%
    mutate(beta_prod = coeff*(-precip_diff)) %>%
    group_by(STATENAME, year, month) %>%
    summarize(change_factor = exp(sum(beta_prod)),
              totaldmg_adj2017 = unique(totaldmg_adj2017),
              cf_dmg = totaldmg_adj2017*change_factor) %>%
    ungroup() %>%
    summarize(obs_damage = sum(totaldmg_adj2017, na.rm = TRUE),
              cf_dmg = sum(cf_dmg, na.rm = TRUE)) %>%
    mutate(diff = obs_damage - cf_dmg,
           model = "MO")

## -----------------------------------------------------------------------------
## SI FIG S3C: PLOT COUNTERFACTUAL DAMAGE COMPARISON
cf_df <- full_join(region_cf_dmg, basin_cf_dmg) %>%
    gather(key, value, cf_dmg, diff) %>%
    mutate(key = factor(key, levels = c("cf_dmg", "diff")),
           model = factor(model, levels = c("region", "MO")))

fig_S3C <- ggplot(cf_df, aes(model, value, fill = key)) +
    geom_col(width = 0.5, alpha = 0.3) +
    scale_y_continuous(name = "Billions (2017$)",
                      breaks = seq(0, 12e9, 3e9),
                      labels = seq(0, 12, 3)) +
    scale_x_discrete(labels = c("using main\nregional model",
                                "accounting for\nupstream basin\nprecipitation"),
                     expand = expansion(add = 0.4)) +
    scale_fill_manual(values = c( "slategray4", "forestgreen"),
                      labels = c("Total observed\ndamages\n",
                                 "Estimated portion\nfrom precipitation\nchange")) +
    geom_text(data = subset(cf_df, key == "diff"), size = 2.5,
              aes(label = paste0("$", round(value/1e9, 1), "B"),
                  y = value + 6e8)) +
    geom_text(size = 2.5, aes(label = "$11.9B",
                  y = 11.9e9 + 6e8)) +
    theme_SI +
    ggtitle("Cumulative Missouri flood damages (1988-2017)") +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.key  = element_rect(size = 3, color = "white"),
          legend.key.size = unit(10, "pt"))


## -----------------------------------------------------------------------------
## COMBINE FIGURES
pdf("./fig_S3.pdf", 6.5, 4.2)
print(ggarrange(
    ggarrange(NULL, fig_S3A, NULL,
              ncol = 3, nrow = 1, widths = c(0.4, 1, 0.4), vjust = 2,
              labels = c("", "A", ""), font.label = list(size = 11, face = "plain")),
    NULL,
    ggarrange(fig_S3B, NULL, fig_S3C,
              ncol = 3, nrow = 1, widths = c(1, .1, 1.1),
              labels = c("B", "", "C"), font.label = list(size = 11, face = "plain")),
    ncol = 1, nrow = 3, heights = c(2, .2, 2)))
dev.off()
