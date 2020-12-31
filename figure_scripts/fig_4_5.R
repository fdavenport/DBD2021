library(tidyverse)
library(raster)
library(sf)
library(ggpubr)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
## FUNCTIONS:
## -----------------------------------------------------------------------------

state_sf <- st_read(state_boundary_file) %>%
    mutate(STATENAME = NAME)

## PLOT RATIO MAP
calc_ratio_map <- function(data, base_thr, years, forcings, quants,
                           pr_varname = "monthly_precip_NORM"){

    ratio_df <- data %>%
        subset(year %in% years & FORCING %in% forcings) %>%
        crossing(quant = quants) %>%
        left_join(base_thr, by = c("centroid_coord", "quant", "modelrun")) %>%
        group_by(centroid_coord, quant, modelrun) %>%
        ## calculate ratios for paired simulations
        summarize(exceed = length(which(get(pr_varname) >= q_thr)),
                  n = length(get(pr_varname)),
                  expected = n*(1-unique(quant)),
                  ri_ratio = exceed/expected) %>%
        ## take ensemble ratio
        ## NOTE: taking the mean of ratios is not the same
        ## sum all actual vs. expected exceedences and then calculate ratio
        group_by(centroid_coord, quant) %>%
        summarize(exceed = sum(exceed),
                  expected = sum(expected),
                  n_above_1 = length(which(ri_ratio > 1)),
                  ri_ratio = exceed/expected) %>%
        dplyr::select(-exceed, -expected)

    ratio_breaks <- c(0.8, 0.9, 1, 1.1, 1.25, 1.5, 2)
    ratio_sf <- raster_df_to_sf(
        ratio_df %>%
        dplyr::select(-n_above_1) %>%
        mutate(quant = paste0("q", quant)) %>%
        group_by(quant, centroid_coord) %>%
        spread(quant, ri_ratio)) %>%
        gather(quant, ri_ratio, contains("q0")) %>%
        mutate(ratio_bin = cut(ri_ratio,
                               breaks = ratio_breaks))

    dot_df <- ratio_df %>%
        separate(centroid_coord, into = c("x", "y"), sep = "_") %>%
        st_as_sf(coords = c("x", "y")) %>%
        st_set_crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
        mutate(p_agree = ifelse(ri_ratio > 1, n_above_1/length(unique(data$modelrun)),
                                1-n_above_1/length(unique(data$modelrun))),
               symbol = ifelse(p_agree > 0.66, NA, 16),
               quant = paste0("q", quant))

    ratio_map <- add_map_formatting(
        ggplot(state_sf) +
        geom_sf(data = ratio_sf, col = "transparent", aes(fill = ratio_bin)) +
        geom_sf(color = "gray36", fill = alpha("white", 0), size = 0.2) +
        geom_sf(data = dot_df, size = 0.35, color = "gray25", aes(shape = symbol)) +
        facet_grid(quant~., switch = "y", labeller=quant_labeller),
        scale_name = "probability ratio",
        showlim = TRUE,
        val = brewer.pal(8, "BrBG")[3:8],
        drp = FALSE) +
        scale_shape_identity()

    return(ratio_map)
}

## -----------------------------------------------------------------------------
## ANALYSIS:
## -----------------------------------------------------------------------------

## find subset of IPCC model realizations that have match in climdex
climdex_files <- list.files(rx5day_cmip_path)
climdex_runs <- do.call(rbind, strsplit(climdex_files, "_"))[,1:3]
colnames(climdex_runs) <- c("MODEL", "FORCING", "RUN")

## for historical, use all historical+rcp8.5 matches
historical_runs <- as.data.frame(climdex_runs) %>%
    mutate(t=1) %>%
    spread(FORCING, t, fill = 0) %>%
    dplyr::select(-rcp26) %>%
    mutate(sum = historical + rcp85) %>%
    subset(sum == 2) %>%
    unite(modelrun, MODEL, RUN) %>%
    subset(modelrun %in% ipcc_runs) %>%
    ## remove simulations with rx5day starting in 1950
    subset(!modelrun %in% c("BNU-ESM_r1i1p1", "FGOALS-g2_r1i1p1")) %>%
    dplyr::pull(modelrun)

print("Historical simulations:")
print(historical_runs)

## for future, use rcp8.5 + rcp2.6 matches
future_runs <- as.data.frame(climdex_runs) %>%
    mutate(t=1) %>%
    spread(FORCING, t, fill = 0) %>%
    mutate(sum = historical + rcp85 + rcp26) %>%
    subset(sum == 3)  %>%
    unite(modelrun, MODEL, RUN) %>%
    subset(modelrun %in% ipcc_runs) %>%
    subset(!modelrun %in% c("BNU-ESM_r1i1p1", "FGOALS-g2_r1i1p1")) %>%
    dplyr::pull(modelrun)

print("Future simulations:")
print(future_runs)

## -----------------------------------------------------------------------------
## READ MONTHLY CMIP5 DATA:
runs <- separate(data.frame(runs = historical_runs), runs, into = c("model", "run"), sep = "_")
file_patterns <- c(paste(paste0(runs$model, "_historical_", runs$run), collapse = "|"),
                   paste(paste0(runs$model, "_rcp85_", runs$run), collapse = "|"),
                   paste(paste0(runs$model, "_rcp26_", runs$run), collapse = "|"))

files <- list.files(path = monthly_cmip_path,
                    pattern = paste(file_patterns, collapse = "|"),
                          full.names = TRUE)

monthly_data <- bind_rows(lapply(files, readRDS)) %>%
    mutate(year = as.numeric(year)) %>%
    unite(modelrun, MODEL, RUN)

## STANDARDIZE PRECIP
model_stats <- monthly_data %>%
    subset(year %in% base_period & FORCING == "historical") %>%
    group_by(modelrun, centroid_coord) %>%
    summarize(mean_precip = mean(monthly_precip),
              sd_precip = sd(monthly_precip))

monthly_data_std <- monthly_data %>%
    left_join(model_stats) %>%
    mutate(monthly_precip_NORM = (monthly_precip - mean_precip)/sd_precip)

rm(monthly_data)

centroid_coords <- unique(monthly_data_std$centroid_coord)

## -----------------------------------------------------------------------------
## READY RX5DAY CMIP5 DATA (CLIMDEX)
rx5day_files <- list.files(path = rx5day_cmip_path,
                    pattern = paste(file_patterns, collapse = "|"),
                    full.names = TRUE)

rx5day_data <- bind_rows(lapply(rx5day_files, readRDS)) %>%
    mutate(year = as.numeric(year),
           centroid_coord = paste0(lon, "_", lat)) %>%
    unite(modelrun, MODEL, RUN) %>%
    ## there are some NAs in the HadGEM2 runs before 1860
    subset(!(modelrun %in% c("HadGEM2-CC_r1i1p1", "HadGEM2-ES_r2i1p1") &
             is.na(rx5day))) %>%
    subset(centroid_coord %in% centroid_coords)

## STANDARDIZE PRECIP
rx5day_stats <- rx5day_data %>%
    subset(year %in% base_period & FORCING == "historical") %>%
    group_by(modelrun, centroid_coord) %>%
    summarize(mean_precip = mean(rx5day),
              sd_precip = sd(rx5day))

rx5day_data_std <- rx5day_data %>%
    left_join(rx5day_stats) %>%
    mutate(rx5day_NORM = (rx5day - mean_precip)/sd_precip)

rm(rx5day_data)

## -----------------------------------------------------------------------------
## CALCULATE BASELINE THRESHOLDS

quants <- c(0.5, 0.75, 0.95, 0.99)
ref_period <- 1860:1920

base_thr_monthly <- monthly_data_std %>%
    subset(year %in% ref_period & FORCING == "historical") %>%
    crossing(quant = quants) %>%
    group_by(centroid_coord, quant, modelrun) %>%
    summarize(q_thr = quantile(monthly_precip_NORM, unique(quant)))

base_thr_rx5day <- rx5day_data_std %>%
    subset(year %in% ref_period & FORCING == "historical") %>%
    crossing(quant = quants) %>%
    group_by(centroid_coord, quant, modelrun) %>%
    summarize(q_thr = quantile(rx5day_NORM, unique(quant)))

monthly_historical <- calc_ratio_map(monthly_data_std, base_thr_monthly,
                    years = 1988:2017,
                    forcings = c("historical", "rcp85"),
                    pr_varname = "monthly_precip_NORM",
                    quants = quants) +
    ggtitle("CMIP5 Historical + RCP8.5\nmonthly precipitation") +
    theme(panel.grid.major=element_line(colour="transparent"),
          strip.text = element_text(angle = 270, colour="black", size = 7))

rx5day_historical <- calc_ratio_map(rx5day_data_std, base_thr_rx5day,
                    years = 1988:2017,
                    forcings = c("historical", "rcp85"),
                    pr_varname = "rx5day_NORM",
                    quants = quants) +
    ggtitle("CMIP5 Historical + RCP8.5\nmonthly max 5-day precipitation")

pdf("./fig_4.pdf", 3.5, 4.45)
print(annotate_figure(ggarrange(monthly_historical,
                                rx5day_historical,
                                ncol = 2, nrow = 1,
                      labels = c("A", "B"),
                      font.label = list(size = 11, face = "plain")),
                      top = text_grob("Change in probability of exceeding\n1860-1920 precipitation thresholds during 1988-2017",
                      size = 7, face = "bold")))
dev.off()


## -----------------------------------------------------------------------------
## FUTURE MAPS
## -----------------------------------------------------------------------------

quants_fut <- c(0.5, 0.95, 0.99)

## subset data to future runs
fut_monthly_std <- monthly_data_std %>%
    subset(modelrun %in% future_runs)

fut_rx5day_std <- rx5day_data_std %>%
    subset(modelrun %in% future_runs)

fut1_period <- 2046:2065
fut2_period <- 2081:2100
calc_change_map <- function(data, years, forcings, quants,
                           pr_varname = "monthly_precip_NORM"){

    base_thr <- data %>%
        subset(year %in% 1988:2017 & FORCING %in% forcings) %>%
        crossing(quant = quants) %>%
        group_by(centroid_coord, quant, modelrun) %>%
        summarize(q_thr_base = quantile(get(pr_varname), unique(quant)))

    change_df <- data %>%
        subset(year %in% years & FORCING %in% forcings) %>%
        crossing(quant = quants) %>%
        group_by(centroid_coord, quant, modelrun) %>%
        summarize(q_thr = quantile(get(pr_varname), unique(quant))) %>%
        left_join(base_thr, by = c("centroid_coord", "quant", "modelrun")) %>%
        mutate(change = q_thr - q_thr_base) %>%
        group_by(centroid_coord, quant) %>%
        summarize(mean_change = mean(change),
                  n_change_pos = length(which(change > 0))) %>%
        mutate(quant = paste0("q", quant))

    change_breaks <- c(-1.5, -1, -0.75, -0.5, -.25, 0, .25, 0.5, 0.75,  1, 1.5)

    change_sf <- raster_df_to_sf(
        change_df %>%
        dplyr::select(-n_change_pos) %>%
        group_by(quant, centroid_coord) %>%
        spread(quant, mean_change)) %>%
        gather(quant, mean_change, contains("q0")) %>%
        mutate(change_bin = cut(mean_change,
                               breaks = change_breaks))

    dot_df <- change_df %>%
        separate(centroid_coord, into = c("x", "y"), sep = "_") %>%
        st_as_sf(coords = c("x", "y")) %>%
        st_set_crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
        mutate(p_agree = ifelse(mean_change > 0, n_change_pos/length(unique(data$modelrun)),
                                1-n_change_pos/length(unique(data$modelrun))),
               symbol = ifelse(p_agree > 0.66, NA, 16))

    change_map <- add_map_formatting(
        ggplot(state_sf) +
        geom_sf(data = change_sf, col = "transparent", aes(fill = change_bin)) +
        geom_sf(color = "gray36", fill = alpha("white", 0), size = 0.2) +
        geom_sf(data = dot_df, size = 0.35, color = "gray25", aes(shape = symbol)) +
        facet_grid(quant~., switch = "y", labeller=quant_labeller),
        scale_name = "change (std. dev.) ",
        showlim = FALSE,
        lab = c("-1", "", "-0.5", "", "0", "", "0.5", "", "1"),
        val = c(brewer.pal(10, "BrBG")[1:5], brewer.pal(10, "RdBu")[6:10]),
        drp = FALSE) +
        scale_shape_identity()

    return(change_map)
}

rcp85_fut1_monthly <- calc_change_map(fut_monthly_std,
                                      years = fut1_period,
                                      forcings = c("historical", "rcp85"),
                                      quants = quants_fut) +
    ggtitle("RCP8.5")

rcp85_fut1_rx5day <- calc_change_map(fut_rx5day_std,
                                      years = fut1_period,
                                      forcings = c("historical", "rcp85"),
                                     quants = quants_fut,
                                     pr_varname = "rx5day_NORM") +
    ggtitle("RCP8.5")

rcp85_fut2_monthly <- calc_change_map(fut_monthly_std,
                                      years = fut2_period,
                                      forcings = c("historical", "rcp85"),
                                      quants = quants_fut) +
    ggtitle("RCP8.5")

rcp85_fut2_rx5day <- calc_change_map(fut_rx5day_std,
                                      years = fut2_period,
                                      forcings = c("historical", "rcp85"),
                                     quants = quants_fut,
                                     pr_varname = "rx5day_NORM") +
    ggtitle("RCP8.5")

## -----------------------------------------------------------------------------
## RCP 2.6

rcp26_fut1_monthly <- calc_change_map(fut_monthly_std,
                                      years = fut1_period,
                                      forcings = c("historical", "rcp26"),
                                      quants = quants_fut) +
    ggtitle("RCP2.6")  +
    theme(panel.grid.major=element_line(colour="transparent"),
          strip.text = element_text(angle = 270, colour="black", size = 7))

rcp26_fut1_rx5day <- calc_change_map(fut_rx5day_std,
                                      years = fut1_period,
                                      forcings = c("historical", "rcp26"),
                                     quants = quants_fut,
                                     pr_varname = "rx5day_NORM") +
    ggtitle("RCP2.6")

rcp26_fut2_monthly <- calc_change_map(fut_monthly_std,
                                      years = fut2_period,
                                      forcings = c("historical", "rcp26"),
                                      quants = quants_fut) +
    ggtitle("RCP2.6") +
    theme(panel.grid.major=element_line(colour="transparent"),
          strip.text = element_text(angle = 270, colour="black", size = 7))

rcp26_fut2_rx5day <- calc_change_map(fut_rx5day_std,
                                      years = fut2_period,
                                      forcings = c("historical", "rcp26"),
                                     quants = quants_fut,
                                     pr_varname = "rx5day_NORM") +
    ggtitle("RCP2.6")


pdf("./fig_5.pdf", 7.2, 3.5)
print(ggarrange(annotate_figure(ggarrange(rcp26_fut2_monthly,
                                          rcp85_fut2_monthly,
                                          ncol = 2, nrow = 1,
                                          labels = c("A", "B"),
                                          font.label = list(size = 11,
                                                            face = "plain")),
                                top = text_grob("Projected change in monthly precipitation by 2081-2100",
                                                size = 8, face = "bold")),
                annotate_figure(ggarrange(rcp26_fut2_rx5day,
                                          rcp85_fut2_rx5day,
                                          ncol = 2, nrow = 1,
                                          labels = c("C", "D"),
                                          font.label = list(size = 11,
                                                            face = "plain")),
                                top = text_grob("Projected change in monthly max 5-day precipitation by 2081-2100",
                                                size = 8, face = "bold"))))
dev.off()


## -----------------------------------------------------------------------------
## SI FIG S7

pdf("./fig_S7.pdf", 7.2, 3.5)
print(ggarrange(annotate_figure(ggarrange(rcp26_fut1_monthly,
                                          rcp85_fut1_monthly,
                                          ncol = 2, nrow = 1),
                                top = text_grob("Projected change in monthly precipitation by 2046-2065",
                                                size = 8, face = "bold")),
                annotate_figure(ggarrange(rcp26_fut1_rx5day,
                                          rcp85_fut1_rx5day,
                                          ncol = 2, nrow = 1),
                                top = text_grob("Projected change in monthly max 5-day precipitation by 2046-2065",
                                                size = 8, face = "bold"))))
dev.off()


