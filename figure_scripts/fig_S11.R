library(dplyr)
library(lfe)
library(ggplot2)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

data <- readRDS(state_panel_file)

## -----------------------------------------------------------------------------

prob_model <- felm(damage_binary ~ monthly_precip_NORM | STATENAME + STATEYEAR + STATEMONTH, data = data)
slope <- prob_model$coefficients[1,]

prob_boot_results <- readRDS("../processed_data/bootstrapped_models/prob_bootstrap_1.Rds")
prob_boot_coeff <- sapply(prob_boot_results, function(x) x[,"Estimate"])
prob_slope_low <- quantile(prob_boot_coeff, 0.025)
prob_slope_hi <- quantile(prob_boot_coeff, 0.975)

xmean <- mean(data$monthly_precip_NORM)
ymean <- mean(data$damage_binary)

x <- seq(-3, 5.5, 0.1)
prob_df <- data.frame(x, dummy = 1) %>%
    mutate(xdiff = x - xmean,
           y = slope*xdiff + ymean,
           yHi = prob_slope_hi*xdiff + ymean,
           yLow = prob_slope_low*xdiff + ymean)


logit_ml <- readRDS("../processed_data/regression_models/logit_model.Rds")
logit_boot <- read_boot_results("../processed_data/bootstrapped_models/logit_model_bootstrap_1.Rds",
                                n_var = 1,
                                var_codes = c("precip")) %>%
    dplyr::rename(coeff = Estimate) %>%
    dplyr::select(-`Std. Error`)

logit_df <- crossing(x, coeff = logit_boot$coeff) %>%
    mutate(logits = x*coeff,
           prob = exp(logits)/(1 + exp(logits))) %>%
    group_by(x) %>%
    dplyr::summarize(low = quantile(prob, 0.025),
                     high = quantile(prob, 0.975)) %>%
    ungroup() %>%
    mutate(x = as.numeric(as.character(x)),
           logits = x*logit_ml$coefficients,
           prob = exp(logits)/(1 + exp(logits)))
logit_shift <- xmean - logit_df[which.min(abs(logit_df$prob - ymean)),]$x


ribbon_color <- brewer.pal(9, "Oranges")[5]
fig_S11 <- ggplot(prob_df, aes(x = x, y = y)) +
    geom_point(data = data, alpha = 0.01, shape = "|", size = 3,
               aes(x = monthly_precip_NORM, y = damage_binary)) +
    geom_ribbon(alpha = 0.6, fill = ribbon_color, aes(ymin = yLow, ymax = yHi)) +
    geom_line() +
    geom_ribbon(data = logit_df, alpha = 0.6, fill = "lightskyblue",
                aes(x = x + logit_shift, y=NULL, ymin = low, ymax = high)) +
    geom_line(data = logit_df, aes(x = x+logit_shift, y = prob), linetype = "dashed",
              size = 0.4) +
    scale_y_continuous(name = "Probability that flood damages occur", expand = c(0,0)) +
    scale_x_continuous(name = "Precipitation anomaly (s.d.)", expand = c(0,0)) +
    coord_cartesian(ylim = c(0,1), xlim = c(-2.5, 5)) +
    theme_pub

pdf("./fig_S11.pdf", 3, 2.5)
print(fig_S11)
dev.off()
