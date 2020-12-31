library(lfe)
library(stargazer)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

data <- readRDS(state_panel_file) %>%
    ungroup()

## -----------------------------------------------------------------------------
## SI FIG S10A: POOLED MODEL FIXED EFFECTS SUMMARY TABLE
linear_model <- readRDS("../processed_data/regression_models/linear.Rds")

pooled_nofe <- felm(damage_value ~ monthly_precip_NORM,
                    data = data)

pooled_sy_fe <- felm(damage_value ~ monthly_precip_NORM | STATEYEAR,
                         data = data)

pooled_s_y_fe <- felm(damage_value ~ monthly_precip_NORM |
                          STATENAME + year,
                      data = data)

pooled_s_y_sm_fe <- felm(damage_value ~ monthly_precip_NORM |
                             STATENAME + year + STATEMONTH,
                         data = data)

syfe <- c("State-year FE", "Yes", "No", "Yes", "No", "No")
sfe <- c("State FE", "No", "No", "No", "Yes", "Yes")
yfe <- c("Year FE", "No", "No", "No", "Yes", "Yes")
smfe <- c("State-month FE", "Yes", "No", "No", "No", "Yes")

pool_list <- list(linear_model, pooled_nofe, pooled_sy_fe,
                  pooled_s_y_fe, pooled_s_y_sm_fe)

pool_rsq_full <- c("R$^2$ (full)",
                   sapply(pool_list, function(x) round(summary(x)$r.squared,3)))
pool_rsq_proj <- c("R$^2$ (proj.)",
                   sapply(pool_list, function(x) round(summary(x)$P.r.squared,3)))
pool_dof <- c("DOF", sapply(pool_list, function(x) summary(x)$df[1]))

stargazer(pool_list,
          title = "Fixed-effects robustness table",
          out = "./tables/pooled_fe_robust.tex",
          column.labels = c("main model"),
          dep.var.labels.include = FALSE,
          dep.var.caption = "",
          covariate.labels = "monthly precip (s.d.)",
          column.sep.width = "2pt",
          keep = "monthly_precip_NORM",
          keep.stat = c("n"),
          omit.stat = "ser",
          df = TRUE,
          font.size = "scriptsize",
          no.space = TRUE,
          add.lines = list(syfe, sfe, yfe, smfe,
                           pool_rsq_full,
                           pool_rsq_proj, pool_dof))

## -----------------------------------------------------------------------------
## SI FIG S10B: POISSON VS OLS

linear_boot_coeff <- read_boot_results("../processed_data/bootstrapped_models/linear_bootstrap_1.Rds",
                                       n_var = 1,
                                       var_codes = "precip") %>%
    dplyr::rename(coeff = Estimate) %>%
    dplyr::select(-`Std. Error`)
print("log-linear model coefficient 95% CI:")
print(quantile(linear_boot_coeff$coeff, c(0.025, 0.975)))

poisson_model <- readRDS("../processed_data/regression_models/poisson_model.Rds")
poisson_boot <- read_boot_results("../processed_data/bootstrapped_models/poisson_model_bootstrap_1.Rds",
                                n_var = 1,
                                var_codes = c("precip")) %>%
    dplyr::rename(coeff = Estimate) %>%
    dplyr::select(-`Std. Error`)
print("Poisson model coefficient 95% CI:")
print(quantile(poisson_boot$coeff, c(0.025, 0.975)))

## -----------------------------------------------------------------------------
## SI FIG S10C: REGIONAL MODEL SUMMARY TABLE
region_model <- readRDS("../processed_data/regression_models/region_model.Rds")

table_var_names <- paste0(gsub("monthly_precip_NORM:", "monthly precip (",
                               gsub("region", "",
                                    rownames(summary(region_model)$coefficients))), ")")

rsq_full <- c("R$^2$ (full)",
              round(summary(region_model)$r.squared,3))
rsq_proj <- c("R$^2$ (proj.)",
              round(summary(region_model)$P.r.squared,3))
region_dof <- c("DOF",
              summary(region_model)$df[1])
r_syfe <- c("State-year FE", "Yes")
r_smfe <- c("State-month FE", "Yes")

stargazer(list(region_model),
          title = "Regional model regression table",
          out = "./tables/main_models.tex",
          column.labels = c("main regional model"),
          dep.var.labels.include = FALSE,
          dep.var.caption = "",
          covariate.labels = table_var_names,
          column.sep.width = "2pt",
          keep.stat = c("n"),
          omit.stat = "ser",
                    df = TRUE,
          font.size = "scriptsize",
          no.space = TRUE,
          add.lines = list(r_syfe, r_smfe, rsq_full, rsq_proj, region_dof))


## -----------------------------------------------------------------------------
## SI FIG S10D: DIAGNOSTIC PLOTS
data_demean <- demeanlist(dplyr::select(data, monthly_precip_NORM, damage_value),
                          list(factor(data$STATEYEAR), factor(data$STATEMONTH))) %>%
    mutate(region = data$region)

main_lm <- lm(damage_value ~ monthly_precip_NORM + monthly_precip_NORM:region,
              data = data_demean)

pdf("./fig_S10.pdf", 7, 1.8)
par(mfrow = c(1,4), mgp = c(2, 1, 0), mar = c(3, 3, 2, 1))
plot(main_lm, sub.caption = "Regression model diagnostic plots",
     cook.levels = NULL,
     cex = 0.2, pch = 20, col = alpha("black", 0.5),
     id.n = 0,
     cex.caption = 0.7,
     cex.axis = 0.8,
     las = 1)
dev.off()
