library(tidyverse); library(RColorBrewer);
library(raster); library(sf);
library(ggpubr); library(quantreg);
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

precip_data <- readRDS(precip_std_file) %>%
    subset(year >= precip_start & year <= precip_end) %>%
    mutate(year_count = year - precip_start + 1)

## read in regional coefficients
region_model <- readRDS("../processed_data/regression_models/region_model.Rds")
region_coeff <- data.frame("coeff" = summary(region_model)$coefficients[,"Estimate"]) %>%
    rownames_to_column(var = "vars") %>%
    separate(vars, sep = ":", into = c("var1", "region")) %>%
    dplyr::select(-var1) %>%
    mutate(region = gsub("region", "", region))

quants <- c(0.5, 0.75, 0.95)
## calculate state p trends
state_p_trends <- precip_data %>%
    crossing(quant = quants) %>%
    group_by(STATENAME, quant) %>%
    summarize(rq_trend_mag = rq(monthly_precip_NORM ~ year,
                                tau = unique(quant))$coefficients[2]) %>%
    ungroup() %>%
    mutate(quant = paste0("q", quant),
           precip_diff = rq_trend_mag*length(precip_start:precip_end)) %>%
    mutate(region = findNCEIregion(STATENAME)) %>%
    merge(region_coeff, by = "region")

boot_trends <- readRDS(random_trend_file) %>%
    subset(range == precip_range)
trend_p_values <- boot_trends %>%
    gather(quant, trend, starts_with("q0")) %>%
    merge(state_p_trends) %>%
    group_by(STATENAME, quant) %>%
    summarize(p_obs = ecdf(trend)(unique(rq_trend_mag)),
              p_0 = ecdf(trend)(0)) %>%
    ungroup() %>%
    mutate(sig = ifelse(p_obs <= 0.025 | p_obs >= 0.975, "sig1",
                 ifelse(p_obs <= 0.05 | p_obs >= 0.95, "sig2",
                        "none")))

statebounds <- shapefile(state_boundary_file)
statebounds$STATENAME <- statebounds$NAME

state_sf <- st_as_sf(statebounds) %>%
    merge(state_p_trends, by = "STATENAME") %>%
    merge(trend_p_values, by = c("STATENAME", "quant")) %>%
    mutate(projected_damage = exp(precip_diff*coeff),
           precip_diff_bin = cut(precip_diff, breaks = c(-Inf, seq(-1, 1, .25), Inf)),
           proj_dmg_bin = cut(projected_damage,
                              breaks = c(-Inf, .3, .5, .8, 1, 1.25, 2, 3, Inf)))

state_sig1 <- subset(state_sf, sig == "sig1")
state_sig2 <- subset(state_sf, sig == "sig2")
sig_col <- c("sig1" = "black", "sig2" = "gray28")
sig_linetype <- c("sig1" = "solid", "sig2" = "dashed")

print("Mean expected damage ratio at each quantile:")
print(as.data.frame(state_sf) %>% group_by(quant) %>% summarize(mean_dmg_ratio = mean(projected_damage, na.rm = TRUE)))
print("Highest ratios of projected damage at the 0.95 quantile:")
print(head(as.data.frame(state_sf) %>% subset(quant == "q0.95") %>%
           arrange(desc(projected_damage)) %>%
           dplyr::select(quant, STATENAME, projected_damage)))

## -----------------------------------------------------------------------------
## SI FIG S4A: MAPS OF PRECIP TRENDS

fig_S4A <- add_map_formatting(
    ggplot(state_sf, aes(fill = precip_diff_bin)) +
    geom_sf(color = "gray60", alpha = 0.8, size = 0.1) +
    geom_sf(data = state_sig2, fill = alpha("white", 0), size = 0.3,
            aes(linetype = "sig2", col = "sig2")) +
    geom_sf(data = state_sig1, fill = alpha("white", 0), size = 0.3,
            aes(linetype = "sig1", col = "sig1")) +
    scale_linetype_manual(name = "sig", values = sig_linetype,
                          labels = c("p<0.05", "p<0.1")) +
    scale_color_manual(name = "sig", values = sig_col,
                     labels = c("p<0.05", "p<0.1")) +
    facet_grid(quant~., switch = "y", labeller=quant_labeller),
    scale_name = "\nstandard deviations",
    val = brewer.pal(10, "BrBG"),
    lab = c("-1", "", "-0.5", "", "0", "", "0.5", "", "1"),
    drp = FALSE) +
    ggtitle("Precipitation change 1928-2017\n") +
    theme(panel.grid.major=element_line(colour="transparent"),
          strip.text = element_text(angle = 270, size = 7, color = "black", vjust = -.25,
                                    margin = margin(r = 3, unit = "pt"))) +
    guides(linetype = "none", color = "none")


sig_legend <- get_legend(fig_S4A +
                 theme(legend.text = element_text(size = 6),
                       legend.title = element_blank()) +
                 guides(fill = "none",
                        linetype = guide_legend(keyheight = unit(1, "pt"),
                                                keywidth = unit(4, "pt"),
                                                direction = "vertical"),
                        color = guide_legend(keyheight = unit(1, "pt"),
                                             keywidth = unit(4, "pt"),
                                             direction = "vertical",
                                             override.aes = list(size = 2))))


## -----------------------------------------------------------------------------
## SI FIG S4B : MAPS OF PROJECTED DAMAGE CHANGE

fig_S4B <- add_map_formatting(
    ggplot(state_sf, aes(fill = proj_dmg_bin)) +
    geom_sf(color = "gray60", size = 0.1) +
    facet_grid(quant~., switch = "y", labeller=quant_labeller),
    scale_name= "ratio of flood damages",
    val = brewer.pal(8, "RdBu"),
    lab = c("0.3", "0.5", "0.8", "1", "1.25", "2", "3"),
    drp = FALSE) +
    ggtitle("Estimated change in flood damages\ndue to precipitation change")

## -----------------------------------------------------------------------------
## COMBINED FIGURE S4

pdf("./fig_S4.pdf", 3.5, 3.5)
print(ggarrange(fig_S4A +
                annotation_custom(sig_legend, xmin = -1700000, xmax = -1600000,
                                 ymin = -2100000, -2000000),
                fig_S4B,
                ncol = 2, nrow = 1,
                hjust = 0,
                labels = c("A", "B"),
                font.label = list(size = 11, face = "plain")))
dev.off()


