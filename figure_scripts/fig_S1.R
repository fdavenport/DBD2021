library(tidyverse)
library(scales)
library(ggpubr)
library(lfe)
source("../analysis_scripts/file_paths.R")
source("../analysis_scripts/func.R")
source("../analysis_scripts/parameters.R")
source("./theme_func.R")

## -----------------------------------------------------------------------------
sheldus_data_state <- read.csv(damage_file) %>%
    filter(!(State.Name %in% c("ALASKA", "DISTRICT OF COLUMBIA",
                               "HAWAII", "PUERTO RICO"))) %>%
    mutate(State.Name = str_to_title(State.Name)) %>%
    group_by(State.Name, Year, Month) %>%
    summarize(hazards = paste(unique(Hazard), collapse = ","),
              TotalDmg.ADJ.2017 = sum(PropertyDmg.ADJ.2017., CropDmg.ADJ.2017.),
              TotalDmg = sum(PropertyDmg, CropDmg))

## -----------------------------------------------------------------------------
## SI FIG S1A: COMPARE REPORTS PER YEAR WITH NFIP CLAIMS DATASET

sheldus_data_year <- sheldus_data_state %>%
    group_by(Year) %>%
    summarize(count = length(TotalDmg.ADJ.2017))

fig_S1Aa <- ggplot(sheldus_data_year, aes(Year, count)) +
    geom_hline(yintercept = c(100, 200), linetype = "dashed") +
    geom_col(alpha = 0.5, col = "dodgerblue3", fill = "dodgerblue3") +
    theme_SI+
    scale_x_continuous(breaks = seq(1960, 2010, 10)) +
    scale_y_continuous(sec.axis = sec_axis(~(.)/576, labels = scales::percent,
                                           name = "percent of total\nstate-months")) +
    coord_cartesian(xlim = c(1960, 2017), ylim = c(0, 280)) +
    ggtitle("Frequency of Flooding damage in SHELDUS") +
    ylab("number of\nstate-months") +
    theme(plot.title = element_text(hjust = 0),
          axis.title.x = element_blank())

nfip_data_year <- read.csv(nfip_claims_file) %>%
    filter(!(state %in% c("DC", "PR", "AK", "HI", "", "GU", "AS", "VI"))) %>%
    dplyr::select(state, dateofloss, amountpaidonbuildingclaim) %>%
    separate(dateofloss, sep = "-", into = c("Year", "Month", "day")) %>%
    mutate(State.Name = state.name[match(state, state.abb)],
           Year = as.numeric(Year),
           Month = as.numeric(Month)) %>%
    group_by(State.Name, Year, Month) %>%
    summarize(number_of_claims = length(amountpaidonbuildingclaim),
              amount_paid = sum(amountpaidonbuildingclaim)) %>%
    full_join(sheldus_data_state) %>%
    filter(Year <= 2017 & Year >= 1978)  %>%
    mutate(SHELDUS_status = ifelse(TotalDmg.ADJ.2017 == 0 | is.na(TotalDmg.ADJ.2017),"missing", "included")) %>%
    group_by(Year) %>%
    filter(!is.na(number_of_claims)) %>%
    summarize(number_of_claims = sum(number_of_claims))

fig_S1Ab <- ggplot(nfip_data_year, aes(Year, number_of_claims/1000)) +
    geom_hline(yintercept = c(100, 200), linetype = "dashed") +
    geom_col(alpha = 0.5, col = "mediumpurple4", fill = "mediumpurple4") +
    theme_SI +
    scale_x_continuous(breaks = seq(1960, 2010, 10)) +
    coord_cartesian(xlim = c(1960, 2017), ylim = c(0, 280)) +
    ggtitle("Frequency of NFIP claims") +
    ylab("NFIP claims\n(thousands)") +
    theme(plot.title = element_text(hjust = 0),
          axis.title.x = element_blank())


pdf("./fig_S1A.pdf", 3.45, 2.2)
print(ggarrange(fig_S1Aa,
               fig_S1Ab,
          ncol = 1, nrow = 2,
          align = "hv"))
dev.off()

## -----------------------------------------------------------------------------
## SI FIG S1C: COMPARE SHELDUS WITH BILLION DOLLAR DISASTERS

bdd <- read.csv(bdd_file)

sheldus_bdd_match <- bdd %>%
    mutate(begin_month = as.Date(paste0(floor(Begin_Date/100), "01"), format = "%Y%m%d"),
           end_month = as.Date(paste0(floor(End_Date/100), "01"), format = "%Y%m%d"),
           date = begin_month) %>%
    group_by(Name) %>%
    complete(., date = seq.Date(begin_month, end_month, by = "month")) %>%
    fill(States) %>%
    separate_rows(States, sep = ", ") %>%
    ungroup() %>%
    mutate(State.Name = state.name[match(States, state.abb)],
           Year = as.numeric(strftime(date, "%Y")),
           Month = as.numeric(strftime(date, "%m"))) %>%
    dplyr::select(Name, Year, Month, State.Name) %>%
    filter(Year <= 2017) %>%
    left_join(sheldus_data_state, by = c("Year", "Month", "State.Name")) %>%
    group_by(Name) %>%
    summarize(SHELDUS_dmg = sum(TotalDmg.ADJ.2017, na.rm = TRUE)) %>%
    right_join(bdd) %>%
    subset(SHELDUS_Flooding == "Yes")


fig_S1C <- ggplot(sheldus_bdd_match, aes(Total_Cost_Millions/1e3, SHELDUS_dmg/1e9)) +
    geom_point(size = 1.5, aes(col = Disaster, shape = Disaster)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_smooth(method = "lm", col = "dodgerblue3", fill = "cadetblue2") +
    scale_color_manual(name = "Type", values = c("black", "darkgreen")) +
    scale_shape_manual(name = "Type", values = c(16, 17)) +
    scale_y_log10(name = "SHELDUS damage\n(Billions US2017$)") +
    scale_x_log10(name = "Billion Dollar Disaster Estimate (Billions US2019$)") +
    theme_SI+
    theme(legend.position = c(.8, .15),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "black", size = .8),
          legend.margin = margin(0, 2, 0, 0, unit = "pt"),
          plot.title = element_text(size = 9, hjust = 1)) +
    ggtitle("SHELDUS vs. NCEI Billion Dollar Disaster damages")

pdf("./fig_S1C.pdf", 3.3, 2)
print(fig_S1C)
dev.off()

## -----------------------------------------------------------------------------
## SIMULATION FOR B and D

set.seed(123)
## data parameters
ns = 48  #number of states
ny = 20*12  #number of time periods (years*months)
nn = ns*ny  #total num of obs
b_true = 1.3 #'true' log coeff
n_sim = 200 ## number of simulations

b_dat <- matrix(,n_sim, 5) # matrix to store estimated betas
colnames(b_dat) <- c("b", "b1", "b2", "b3", "b4")

for(i in 1:n_sim){

  ## 'true' y data
  dat <- data.frame(year = rep(1:ny,each=ns), state = rep(1:ns,times=ny),
                    x = rnorm(nn,10,1), e = rnorm(nn,0,2))
  dat$y <- exp(dat$x*b_true+dat$e)
  dat$logy <- ifelse(dat$y>0, log(dat$y), log(1))
  b <- coef(felm(logy ~ x | state + year, data=dat))

  ## y data with 10% of events unreported across all years
  dat$y_missing1 <- dat$y
  dat$y_missing1[sample(nn, (nn)*0.1)] <- 0
  dat$logy_missing1 <- ifelse(dat$y_missing1 >0, log(dat$y_missing1), log(1))
  b1 <- coef(felm(logy_missing1 ~ x | state + year, data=dat))

  ## y data with 20% of events unreported across all years
  dat$y_missing2 <- dat$y
  dat$y_missing2[sample(nn, nn*0.2)] <- 0
  dat$logy_missing2 <- ifelse(dat$y_missing2 >0, log(dat$y_missing2), log(1))
  b2 <- coef(felm(logy_missing2 ~ x | state + year, data=dat))

  ## y data with 40% of events unreported in first half of time period
  dat$y_missing3 <- dat$y
  dat$y_missing3[sample(nn/2, (nn/2)*.4)] <- 0
  dat$logy_missing3 <- ifelse(dat$y_missing3 >0, log(dat$y_missing3), log(1))
  b3 <- coef(felm(logy_missing3 ~ x | state + year, data=dat))

  ## y data with all events randomly underestimated (40% of true value on average)
  dat$y_low <- dat$y*runif(nn, 0, .8)
  dat$logy_low <- ifelse(dat$y_low > 0, log(dat$y_low), log(1))
  b4 <- coef(felm(logy_low ~ x | state + year, data=dat))

  b_dat[i,] <- c(b, b1, b2, b3, b4)
}

fig_S1B <- ggplot(as.data.frame(b_dat)) +
  geom_vline(xintercept = b_true, linetype = "dashed") +
  geom_text(x = b_true+.02, y = 28, label = expression("true"~beta),
            hjust = 'left', fontface = "plain", size = 2) +
  stat_density(aes(b, col = "btrue"), geom = 'line', position = 'identity') +
  stat_density(aes(b1, col = "missing1"), geom = 'line', position = 'identity') +
  stat_density(aes(b2, col = "missing2"), geom = 'line', position = 'identity') +
  scale_color_manual(labels = c("'true' data", "10% of events unreported",
                                "20% of events unreported"),
                     values = c("black", "purple", "red")) +
  scale_x_continuous(name = expression("Estimated effect of precipitation on damages ("~beta~")"),
                     lim = c(0.75, 1.55), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 30)) +
  theme_classic(9) +
  theme(legend.position = c(0.33, 0.8),
        legend.title = element_blank(),
        legend.key.size = unit(12, "pt"),
        plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 8)) +
  ggtitle("Simulation where some damages are missing\n(i.e. unreported)")


pdf("./fig_S1B.pdf", 3.3, 2)
print(fig_S1B)
dev.off()

fig_S1D <- ggplot(as.data.frame(b_dat)) +
  geom_vline(xintercept = b_true, linetype = "dashed") +
  geom_text(x = b_true+.02, y = 28, label = expression("true"~beta),
            hjust = 'left', fontface = "plain", size = 2) +
  stat_density(aes(b, col = "btrue"), geom = 'line', position = 'identity') +
  stat_density(aes(b4, col = "low"), geom = 'line', position = 'identity') +
  scale_color_manual(labels = c("'true' data", "damages randomly\nunderestimated"),
                     values = c("black", "dodgerblue2")) +
  scale_x_continuous(name = expression("Estimated effect of precipitation on damages ("~beta~")"),
                     lim = c(0.75, 1.55), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), lim = c(0,30)) +
  theme_classic(9) +
  theme(legend.position = c(0.25, 0.85),
        legend.title = element_blank(),
        legend.key.size = unit(12, "pt"),
        plot.title = element_text(size = 9),
        axis.title.x = element_text(size = 8)) +
  ggtitle("Simulation where value of damage is randomly\nunderestimated")

pdf("./fig_S1D.pdf", 3.3, 2)
print(fig_S1D)
dev.off()


## -----------------------------------------------------------------------------
## TABLE: DISTRIBUTION OF DAMAGES BY REGION

data <- readRDS(state_panel_file)

sheldus_region <- data %>%
    ungroup() %>%
    mutate(tot = sum(totaldmg_adj2017, na.rm = TRUE),
           tot_events = length(which(totaldmg_adj2017 > 0))) %>%
    group_by(region) %>%
    summarize(n_events = length(which(totaldmg_adj2017 > 0)),
              percent_of_events = percent(n_events/unique(tot_events)),
              damage = dollar(sum(totaldmg_adj2017, na.rm = TRUE),
                              largest_with_cents = 1e+03),
              percent_of_damage = percent(sum(totaldmg_adj2017, na.rm = TRUE)/
                                          unique(tot))) %>%
    mutate(region = as.character(region))

## TABLE: DISTRIBUTION OF DAMAGES BY SEASON
sheldus_season <- data %>%
    ungroup() %>%
    mutate(tot = sum(totaldmg_adj2017, na.rm = TRUE),
           tot_events = length(which(totaldmg_adj2017 > 0))) %>%
    group_by(season) %>%
    summarize(n_events = length(which(totaldmg_adj2017 > 0)),
              percent_of_events = percent(n_events/unique(tot_events)),
              damage = dollar(sum(totaldmg_adj2017, na.rm = TRUE),
                              largest_with_cents = 1e+03),
              percent_of_damage = percent(sum(totaldmg_adj2017, na.rm = TRUE)/
                                          unique(tot))) %>%
    mutate(season = as.character(season))

print(sheldus_region)

print(sheldus_season)

## -----------------------------------------------------------------------------
## TABLE: PROPORTION OF MONTHS WITH REPORTED DAMAGE

print("Proportion of state-months without reported damages:")
print(length(which(data$damage_value == 0))/nrow(data))

prop_dmg <- data %>%
    group_by(STATENAME) %>%
    summarize(nzero = length(which(damage_value == 0)),
              prop_zero = nzero/length(damage_value),
              prop_dmg = round(1-prop_zero, 2)) %>%
    mutate(STATENAME = as.character(STATENAME))

overall_prop <- data %>%
    ungroup() %>%
    summarize(nzero = length(which(damage_value == 0)),
              prop_zero = nzero/length(damage_value),
              prop_dmg = round(1-prop_zero, 2)) %>%
    mutate(STATENAME = "All")

print(bind_rows(overall_prop, prop_dmg)[,c("STATENAME", "prop_dmg")], n = 49)

