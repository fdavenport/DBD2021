library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(quantreg)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")
## -----------------------------------------------------------------------------

precip_df <- readRDS(precip_std_file)

## create ranges for precipitation trend calculations
ppt_ranges <- crossing(start = seq(1898, 1988, 5), end = 2017) %>%
    mutate(diff = end - start + 1)
ppt_trends <- vector('list', nrow(ppt_ranges))

for(i in 1:nrow(ppt_ranges)){

    ppt_trends[[i]] <- precip_df %>%
        subset(year >= ppt_ranges$start[i] & year <= ppt_ranges$end[i]) %>%
        crossing(quant = quants) %>%
        group_by(STATENAME, quant) %>%
        summarize(rq_trend_mag = rq(monthly_precip_NORM ~ year,
                                    tau = unique(quant))$coefficients[2]) %>%
        ungroup() %>%
        mutate(quant = paste0("q", quant),
               start = ppt_ranges$start[i],
               end = ppt_ranges$end[i],
               range = paste0(start, "-", end),
               years = end-start+1,
               diff17 = rq_trend_mag*years,
               diff88 = rq_trend_mag*(years-30))

}

ppt_trends <- bind_rows(ppt_trends)

fig_S5A <- ggplot(ppt_trends, aes(start, rq_trend_mag, group = start)) +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_point(alpha = 0.2, size = 0.7) +
    facet_wrap(~quant, nrow = 1, ncol = 3, labeller=quant_labeller,
               scales = "free") +
    scale_x_continuous(name = "Start of time period for trend calculation") +
    scale_y_continuous(name = "Magnitude of trend (s.d. per year)") +
    ggtitle("Sensitivity of precipitation trends to starting year of trend calculation") +
    theme_SI

## -----------------------------------------------------------------------------
## STATE UNCERTAINTY USING DELETE-D RESAMPLING

ppt_trends_delete_d <- readRDS(paste0("../processed_data/quantile_trends_",
                                      precip_start, "-", precip_end,
                                      ".Rds")) %>%
    mutate(quant = paste0("q", q_mid),
           region = factor(findNCEIregion(STATENAME), levels = region_order),
           diff17 = bin_rq_trend*length(trend_start:damage_end)) %>%
    left_join(state_order)

label_df <- ppt_trends_delete_d %>%
    group_by(abb, quant) %>%
    summarize(high = max(diff17)) %>%
    mutate(offset = ifelse(quant == "q0.5", 0.1,
                       ifelse(quant == "q0.75", 0.15, 0.25))) %>%
    ungroup()

fig_S5B <- ggplot(ppt_trends_delete_d,
                                 aes(abb, diff17)) +
    geom_hline(yintercept = 0, size = 0.3) +
    geom_boxplot(size = 0.3, outlier.size = 0.25, outlier.stroke = 0.25,
                 outlier.color = "grey36", aes(fill = region)) +
    geom_text(data = label_df, size = 2.5, angle = 90,
              aes(x = abb, label = abb, y = high + offset)) +
    facet_wrap(~quant, ncol = 1, nrow = 3, scales = "free",
               labeller = quant_labeller) +
    scale_fill_manual(values = alpha(region_colors, 0.7)) +
    scale_y_continuous(name = paste0("Change in precipitation by 2017 (std. dev.) based on ",
                                     precip_start, "-2017 trend")) +
    theme_SI +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom",
          axis.title.x = element_blank())

## -----------------------------------------------------------------------------

pdf("./fig_S5.pdf", 7.2, 7.5)
print(ggarrange(fig_S5A,
                NULL,
                fig_S5B,
                ncol = 1,
                nrow = 3,
                heights = c(1.3, .1, 3),
                labels = c("A","", "B"),
                font.label = list(size = 11, face = "plain")))
dev.off()
