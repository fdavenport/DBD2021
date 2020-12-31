library(tidyverse)
library(ggpubr)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
## SI FIG S9A: STATE-LEVEL INCOME vs. HOUSING VALUE ESTIMATES
## -----------------------------------------------------------------------------

## INCOME DATA - state-year
state_income <- read.csv(income_file) %>%
    dplyr::rename(STATENAME = GeoName) %>%
    dplyr::select(-GeoFips, -LineCode) %>%
    mutate(income_billion_current = income_millions_current/1000)

home_values_state <- read.csv(home_value_state_file)

hu_state <- read.csv(housing_units_file)

state_exposure_df <- home_values_state %>%
    subset(Adjusted == "current") %>%
    full_join(hu_state) %>%
    full_join(state_income) %>%
    mutate(housing_stock_value_billions = home_value*(hu_est/1e9)) %>%
    na.omit()

home_value_model <- lm(log(housing_stock_value_billions) ~ log(income_billion_current),
                       data = state_exposure_df)

fig_S9A <- ggplot(state_exposure_df, aes(income_billion_current,
                                    housing_stock_value_billions)) +
    geom_point(size = 1, alpha = 0.6, aes(shape = as.factor(year))) +
    annotate("text", x = 10, y = 3000,
              label = bquote(r^2==~.(round(summary(home_value_model)$r.squared,3)))) +
    theme_classic_black() +
    coord_cartesian(xlim = c(4,1700), ylim = c(4,5100)) +
    scale_y_log10(name = "Estimated housing value ($ billions)\nmedian home value * number of housing units") +
    scale_x_log10(name = "Annual state income ($ billions)") +
    ggtitle("Housing value vs. income (state-level)") +
    geom_smooth(method = "lm", se = FALSE, col = "#00AD56") +
    theme_SI +
    theme(legend.title = element_blank(),
          legend.position = c(0.8, 0.2))

## -----------------------------------------------------------------------------
## SI FIG S9B: NATIONAL-LEVEL INCOME vs. FIXED ASSETS
## -----------------------------------------------------------------------------

## net stock of reproducible fixed assets (annual)
natl_fixed_assets <- read.csv(assets_file)

natl_income <- state_income %>% group_by(year) %>%
    dplyr::summarize(income_billions_current = sum(income_millions_current, na.rm = TRUE)/1e3)

natl_exposure_df <- full_join(natl_fixed_assets, natl_income, by = "year") %>%
    subset(year >= damage_start & year <= damage_end)

natl_asset_model <- lm(totalAssets_current ~ income_billions_current, data = natl_exposure_df)

fig_S9B <- ggplot(natl_exposure_df, aes(income_billions_current,
                                                   totalAssets_current)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, col = "#00AD56") +
    annotate("text", x = 19000, y = 63000,
             label = bquote(r^2==~.(round(summary(natl_asset_model)$r.squared,3)))) +
    scale_x_continuous(name = "Annual national income ($ billions)") +
    scale_y_continuous(name = "Net stock of reproducible fixed assets ($ billions)") +
    ggtitle("Fixed assets vs. income (country-level)\n1988-2017") +
    theme_SI

## -----------------------------------------------------------------------------

pdf("./fig_S9.pdf", 3, 5.2)
print(ggarrange(fig_S9A, NULL,
                fig_S9B,
                ncol = 1, nrow = 3,
                labels = c("A","",  "B"),
                vjust = c(1.5, 0, 0.5),
                font.label = list(size = 11, face = "plain"),
                heights = c(2.5, .3, 2.5)))
dev.off()
