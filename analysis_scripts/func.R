library(quantreg)
library(sf)
## -----------------------------------------------------------------------------
## FUNCTIONS

findNCEIregion <- function(st){
    ## assign NCEI region
    ## INPUTS: st = statename in all caps (i.e. "ALABAMA")
    ## RETURNS: region name

    swStates <- c("COLORADO", "UTAH", "NEW MEXICO", "ARIZONA")
    wStates <- c("CALIFORNIA", "NEVADA")
    nwStates <- c("WASHINGTON", "OREGON", "IDAHO")
    nrpStates <- c("WYOMING", "MONTANA", "NORTH DAKOTA", "SOUTH DAKOTA", "NEBRASKA")
    southStates <- c("TEXAS", "KANSAS", "OKLAHOMA", "ARKANSAS", "LOUISIANA", "MISSISSIPPI")
    cenStates <- c("MISSOURI", "ILLINOIS", "INDIANA", "OHIO", "KENTUCKY", "WEST VIRGINIA",
                   "TENNESSEE")
    umStates <- c("MINNESOTA", "WISCONSIN", "IOWA", "MICHIGAN")
    neStates <- c("PENNSYLVANIA", "NEW YORK", "VERMONT", "NEW HAMPSHIRE", "MAINE",
                  "MASSACHUSETTS", "RHODE ISLAND", "CONNECTICUT", "NEW JERSEY",
                  "DELAWARE", "MARYLAND")
    seStates <- c("VIRGINIA", "NORTH CAROLINA", "SOUTH CAROLINA", "GEORGIA", "ALABAMA",
                  "FLORIDA")

    st <- toupper(st)
    ifelse(st %in% swStates, "Southwest",
    ifelse(st %in% wStates, "West",
    ifelse(st %in% nwStates, "Northwest",
    ifelse(st %in% nrpStates, "Northern Rockies Plains",
    ifelse(st %in% southStates, "South",
    ifelse(st %in% cenStates, "Central",
    ifelse(st %in% umStates, "Upper Midwest",
    ifelse(st %in% neStates, "Northeast",
    ifelse(st %in% seStates, "Southeast", "otherNCEI")))))))))
}


## -----------------------------------------------------------------------------
## CONVERT DATA FRAME OF RASTER INFORMATION INTO SF OBJECT
raster_df_to_sf <- function(df){
    new_sf <- st_as_sf(rasterToPolygons(rasterFromXYZ(
        df %>%
        separate(centroid_coord, into = c("x", "y"), sep = "_") %>%
    dplyr::select(x, y, everything()))))
    st_crs(new_sf) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    return(new_sf)
}

## -----------------------------------------------------------------------------

format_coeff <- function(dat, n_vars = 1, var_codes, POLY = FALSE){
    ## Format regression coefficient results
    ## INPUTS:
    ##  dat = regression summary coefficients
    ##  n_vars = numeric, number of vars to separate
    ##  var_codes = character vector of variables to match
    ##  POLY = TRUE if model includes poly formula
    ## RETURNS: formatted data.frame of coefficients

    if(class(dat) == "felm") {
        dat <- as.data.frame(summary(dat)$coefficients) %>%
            rownames_to_column(var = "vars")
    }

    tmpvar <- paste0("var", 1:n_vars)

    formatted_dat <- dat %>%
        separate(vars, sep = ":", into = tmpvar)

    for(i in var_codes){
        formatted_dat[,i] <- sapply(1:nrow(formatted_dat),
                                    function(x) {
                                        col_ind <- which(grepl(i, formatted_dat[x,]))
                                        ifelse(length(col_ind) == 0,
                                               NA,
                                               formatted_dat[x,col_ind])})
        if(i != "precip"){
            formatted_dat[,i] <- gsub(i, "", formatted_dat[,i])
        }
    }
    if(POLY){
        formatted_dat <- formatted_dat %>%
            separate(precip, into = c("junk", "precip", "poly_order"),
                     sep = "[(]|[)]") %>%
            separate(precip, into = "precip", sep = ",") %>%
            dplyr::select(-junk)
    }
    formatted_dat <- formatted_dat %>%
        dplyr::select(-starts_with("var"), -starts_with("t value"), -starts_with("Pr(>"),
                       -starts_with("z value"))
    return(formatted_dat)
}

## -----------------------------------------------------------------------------

read_boot_results <- function(files, n_var, var_codes, POLY = FALSE){
    ## Reads and formats bootstrapped regression results
    ## INPUTS:
    ##  files = vector of filenames to read (written from model-bootstrapping script)
    ##  n_vars = numeric, number of vars to separate
    ##  var_codes = character vector of variables to match
    ##  POLY = TRUE if model include poly formula
    ## RETURNS: formatted data.frame of coefficients

    boot_results <- bind_rows(lapply(unlist(lapply(files, readRDS), recursive = FALSE),
                           function(x) {
                               as.data.frame(x) %>%
                                   rownames_to_column(var = "vars")}), .id = "boot_id") %>%
        format_coeff(., n_var, var_codes, POLY)
}

## -----------------------------------------------------------------------------

mbb <- function(orig_dat, block_size = 1, nboot = 1000){
    ## Moving block bootstrap
    ## Samples timeseries using moving block bootstrap
    ## (allows sampling of multiple variables at once)
    ## INPUTS:
    ##  orig_dat = matrix of data to be sampled (variables in columns, observations in rows)
    ##  block_size = size of block for bootstrapping (default = 1)
    ##  nboot = number of bootstraps(default = 1000)
    ## RETURNS: list of length nboot, where each element is resampled version of orig_dat.
    ## If nboot = 1, result is unlisted first

    if(!is.data.frame(orig_dat)) orig_dat <- data.frame(orig_dat)

    n <- nrow(orig_dat)
    if(n %% block_size != 0) {
        warning("Time series length is not integer multiple of block size.",
                "Resample size will be (n %/% block_size)*block_size")
    }

    new_data <- vector('list', length = nboot)
    for(i in 1:nboot){
        block_smp <- sample(1:(n-block_size+1), ceiling(n/block_size), replace = TRUE)
        smp <- NULL
        for(j in seq_along(block_smp)){
            smp <- c(smp, block_smp[j]:(block_smp[j]+block_size-1))
        }
        smp <- smp[1:n]
        new_data[[i]] <- orig_dat[smp,]
    }
    if(nboot == 1) new_data <- unlist(new_data)

    new_data
}


## -----------------------------------------------------------------------------

mbb_rq_ts <- function(x, y, tau, block_size = 1, nboot = 1000){
    ## calculate quantile regression trends for timeseries resampled using mbb
    ## INPUTS:
    ##  x = time variable for quantile regression; vector
    ##  y = dependent variable; vector of same length as x
    ##  tau = desired quantile, e.g. 0.5 for median
    ##  block_size = size of block for bootstrapping (default = 1)
    ##  nboot = number of bootstraps (default = 1000)
    ## RETURNS: nboot x tau matrix with estimated rq coefficients for each bootstrap sample

    resampled_data <- mbb(data.frame(y = y), block_size, nboot)


    slope_results <- matrix(, nrow = nboot, ncol = length(tau))

    for(i in seq_along(resampled_data)){
        if(length(tau) == 1){
            ## calculate trend on new y values vs. original x
            slope_results[i,] <- rq(resampled_data[[i]] ~ x, tau = tau)$coefficients[2]
        } else {
             slope_results[i,] <- rq(resampled_data[[i]] ~ x, tau = tau)$coefficients[2,]
        }
    }

    colnames(slope_results) <- paste0("q", tau)
    slope_results
}

## -----------------------------------------------------------------------------

delete_d_sample <- function(x, d, n = NULL){
    ## Delete-d years from time-series
    ## INPUTS:
    ##  x = years to sample from
    ##  d = number of years to be deleted
    ##  n = number of iterations, optional (default is to calculate all possible iterations)

    if(is.null(n)) n <- choose(length(x), d)

    include_years <- vector('list', length = n)
    for(i in 1:n){
        include_years[[i]] <- sample(x, length(x)-d, replace = FALSE)
    }
    include_years
}

## -----------------------------------------------------------------------------

calc_bin_breaks <- function(dat, q_breaks, allyears, group_vars = "STATENAME"){
    ## Calculate breaks for each quantile bin
    ## INPUTS:
    ##  dat = data.frame with monthly_precip_NORM, year and grouping variables (other variables OK)
    ##  q_breaks = numeric vector of quantile values separating bins
    ##  allyears = numeric vector of years to calculate breaks values for
    ##  group_vars = character vector of variables to group by (default = "STATENAME")
    ## RETURNS: data.frame with quantile bin breaks in each state in each year

    quantile_breaks <- dat %>%
        crossing(q_break = q_breaks) %>%
        group_by_at(vars(c(group_vars, "q_break"))) %>%
        summarize(break_trend = rq(monthly_precip_NORM ~ year,
                                   tau = unique(q_break))$coefficients[2],
                  break_intercept = rq(monthly_precip_NORM ~ year,
                                       tau = unique(q_break))$coefficients[1]) %>%
        crossing(year = allyears) %>%
        mutate(break_value = year*break_trend + break_intercept)
}

## -----------------------------------------------------------------------------

calc_bin_trends <- function(dat, q_mid, group_vars = "STATENAME"){
    ## Calculate trends for each quantile bin
    ## INPUTS:
    ##  dat = data.frame with monthly_precip_NORM, year and grouping variables (other variables OK)
    ##  q_mid = numeric vector of quantile values in the middle of each bin
    ##  group_vars = character vector of variables to group by (default = "STATENAME")
    ## RETURNS: data.frame with trends for each quantile bin

    quantile_trends <- dat %>%
        crossing(q_mid = q_mid) %>%
        group_by_at(vars(c(group_vars, "q_mid"))) %>%
        summarize(bin_rq_trend = rq(monthly_precip_NORM ~ year,
                                    tau = unique(q_mid))$coefficients[2]) %>%
        ungroup()
}

## -----------------------------------------------------------------------------
detrend_precip <- function(dat, quantile_breaks, quantile_trends, trend_start,
                           group_vars = "STATENAME"){
    ## Detrend precip based on quantile trend
    ## INPUTS:
    ##  dat = data.frame with monthly_precip_NORM, year and grouping variables (other variables OK)
    ##  quantile_breaks = data.frame returned from calc_bin_breaks
    ##  quantile_trends = data.frame returned from calc_bin_trends
    ##  trend_start = numeric starting year of trend calculation period
    ##  group_vars = character vector of variables to group by (default = "STATENAME")
    ## RETURNS: original data.frame plus additional variables precip_diff and monthly_precip_NORM_detrend

    q_mid <- unique(quantile_trends$q_mid)
    ## assign each precip value to bin
    ## calculate counterfactual precip
    detrended_precip <- dat %>%
        merge(quantile_breaks) %>%
        group_by(STATENAME, year, month) %>%
        summarize(q_bin = cut(unique(monthly_precip_NORM), c(-10, unique(break_value), 10),
                              labels = FALSE),
                  monthly_precip_NORM = unique(monthly_precip_NORM)) %>%
        mutate(q_mid = q_mid[q_bin]) %>%
        merge(quantile_trends) %>%
        mutate(precip_diff = bin_rq_trend*(year - trend_start),
               monthly_precip_NORM_detrend = monthly_precip_NORM - precip_diff) %>%
        dplyr::select(-q_mid, -q_bin, -bin_rq_trend)
}

## -----------------------------------------------------------------------------

calc_cf_damage <- function(precip_dat, damage_dat, coeffs, LAG = FALSE, QUAD = FALSE){
    ## Calculate counterfactual damages from detrended precipitation time series
    ## INPUTS:
    ##  precip_dat = (bootstrapped) detrended precipitation data.frame
    ##  damage_dat = data.frame with damages data
    ##  coeffs = (bootstrapped) model coefficients data.frame
    ##  LAG = TRUE if lagged model
    ##  QUAD = TRUE if quadratic model
    ## RETURNS: counterfactual damage timeseries (data.frame)

    if(LAG){
        cf_dat <- precip_dat %>%
            dplyr::select(STATENAME, year, month, precip_diff) %>%
            ## CALCULATE LAGGED PRECIP FOR EACH STATE
            group_by(STATENAME) %>%
            arrange(year, month) %>%
            mutate(monthly_precip_NORM_lag1_diff = lag(precip_diff, 1,
                                                       order_by = STATENAME),
                   monthly_precip_NORM_lag2_diff = lag(precip_diff, 2,
                                                       order_by = STATENAME)) %>%
            ungroup() %>%
            right_join(damage_dat) %>%
            ## SPREAD PRECIP VARIABLES AND CALCULATE CHANGE FACTOR
            dplyr::rename(monthly_precip_NORM_diff = precip_diff) %>%
            gather(precip, detrended_diff, contains("_diff")) %>%
            mutate(precip = gsub("_diff", "", precip)) %>%
            left_join(coeffs) %>%
            mutate(beta_prod = coeff*(-detrended_diff)) %>%
            dplyr::select(-detrended_diff, -coeff) %>%
            spread(precip, beta_prod) %>%
            mutate(change_factor = exp(monthly_precip_NORM +
                                       monthly_precip_NORM_lag1 +
                                       monthly_precip_NORM_lag2),
                   cf_dmg = totaldmg_adj2017*change_factor)
    } else if(QUAD){
        cf_dat <- precip_dat %>%
            dplyr::select(STATENAME, year, month, monthly_precip_NORM, precip_diff) %>%
            right_join(damage_dat) %>%
            left_join(coeffs  %>% pivot_wider(names_from = poly_order, names_prefix = "poly",
                                             values_from = coeff)) %>%
            ## CALCUALTE BETA FROM QUADRATIC
            mutate(beta = poly1 + poly2*monthly_precip_NORM,
                   change_factor = exp(beta*(-precip_diff)),
                   cf_dmg = totaldmg_adj2017*change_factor)
    } else {
        cf_dat <- precip_dat %>%
            dplyr::select(STATENAME, year, month, precip_diff) %>%
            right_join(damage_dat) %>%
            left_join(coeffs) %>%
            mutate(beta_prod = coeff*(-precip_diff),
                   change_factor = exp(beta_prod),
                   cf_dmg = totaldmg_adj2017*change_factor)
    }

    return(cf_dat)
}
