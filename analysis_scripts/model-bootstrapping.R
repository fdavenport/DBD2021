library(dplyr)
library(lfe)
library(fixest)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

## -----------------------------------------------------------------------------
## THIS SCRIPT IS WRITTEN FOR GENERAL BOOTSTRAPPING OF A FIXED EFFECTS LINEAR MODEL
## USING FELM or FEPOIS
## ALL DATA-WRANGLING TAKES PLACE BEFOREHAND

## read in arguments from sbatch script
args <- commandArgs(trailingOnly = TRUE)
ARRAY_ID <- as.numeric(args[1]) ## job number within array - used to set seed and in outfile name
N_BOOT_PER_JOB <- as.numeric(args[2]) ## number of bootstraps for each job
MODEL_NAME <- args[3] ## shorthand name to save output
MODEL_FORM <- args[4] ## string specifying model formula
DATA_FILE <- args[5] ## location of dataframe to use
CLUSTER_VARIABLE <- args[6] ## which variable to sample by (ie. state or county)
if(CLUSTER_VARIABLE == "TIME") BLOCK_SIZE <- as.numeric(args[7])

## OUTPUT FILE -----------------------------------------------------------------
outfile <- paste0("../processed_data/bootstrapped_models/", MODEL_NAME,
                  "_bootstrap_", ARRAY_ID, ".Rds")

## ECHO INPUT VARIABLES --------------------------------------------------------
print("Model Name:")
print(MODEL_NAME)
print("Model Formula:")
print(MODEL_FORM)
print("Bootstrapping by...")
print(CLUSTER_VARIABLE)
if(CLUSTER_VARIABLE == "TIME"){
    print("Block size:")
    print(BLOCK_SIZE)
}
print("Output File:")
print(outfile)
## -----------------------------------------------------------------------------

## bootstrap within each region if model form includes regional interaction
REGION <- grepl("region", MODEL_FORM)
POISSON <- grepl("poisson", MODEL_NAME)
LOGIT <- grepl("logit", MODEL_NAME)

## read in data
data <- readRDS(DATA_FILE)

set.seed(100*ARRAY_ID)

if(CLUSTER_VARIABLE == "TIME"){
    years <- as.character(unique(data$year))
    data_split <- split(data, f = data[,"year"])
    sample_list <- mbb(years, block_size = BLOCK_SIZE, nboot = N_BOOT_PER_JOB)

} else if(REGION){
    ## IF THERE IS REGIONAL INTERACTION TERM, RESAMPLE WITHIN EACH REGION
    data_region <- split(data, f = data[,"region"])
    N_region <- sapply(data_region, function(x) length(unlist(unique(x[,CLUSTER_VARIABLE]))))
    cumN_region <- c(0, cumsum(N_region))

    ## split data by states in "region order" to preserve sampling by region
    data_split <- split(data,
                        f = data[,CLUSTER_VARIABLE])[as.character(state_order$STATENAME)]

    sample_list <- vector('list', length = N_BOOT_PER_JOB)
    for(k in 1:N_BOOT_PER_JOB){
        sample_list[[k]] <- unlist(lapply(seq_along(N_region),
                                          function(z) sample(1:N_region[z], N_region[z],
                                                             replace = TRUE) + cumN_region[z]))

    }
} else {
    ## split into list of dataframes (one for each clustered variable)
    data_split <- split(data, f = data[,CLUSTER_VARIABLE])

    ## number of clusters to sample (i.e. original number of counties or states in dataset)
    N <- length(unlist(unique(data[,CLUSTER_VARIABLE])))

    ## generate matrix of random numbers to resample with replacement
    sample_list <- lapply(1:N_BOOT_PER_JOB, function(x) sample(1:N, N, replace = TRUE))
}

## create matrix to store results
boot_results <- vector('list', length = N_BOOT_PER_JOB)
if(MODEL_NAME %in% c("region_model", "region_5day_model")){
    diff_results <- vector('list', length = N_BOOT_PER_JOB)
}

for(i in 1:N_BOOT_PER_JOB){
    if(i%%100 == 0) print(i)
    current_data <- do.call(bind_rows, data_split[sample_list[[i]]])
    if(POISSON){
        current_model <- fepois(formula(MODEL_FORM), data = current_data)
        boot_results[[i]] <- summary(current_model)$coeftable
    } else if(LOGIT){
        current_model <- femlm(formula(MODEL_FORM), family = "logit",
                               data = current_data)
        boot_results[[i]] <- summary(current_model)$coeftable
    } else {
        current_model <- felm(formula(MODEL_FORM), data = current_data)
        boot_results[[i]] <- summary(current_model)$coefficients
    }


    if(MODEL_NAME == "region_model"){
        linear_model <- felm(damage_value ~ monthly_precip_NORM | STATENAME + STATEYEAR,
                             data = current_data)
        diff_results[[i]] <- current_model$coefficients -
            as.numeric(linear_model$coefficients)
    } else if(MODEL_NAME == "region_5day_model"){
         linear_model <- felm(damage_value ~ precip_5daymax_NORM | STATENAME + STATEYEAR,
                             data = current_data)
        diff_results[[i]] <- current_model$coefficients -
            as.numeric(linear_model$coefficients)
    }

}

saveRDS(boot_results, outfile)

if(MODEL_NAME %in% c("region_model", "region_5day_model")){
    saveRDS(diff_results,
            paste0("../processed_data/bootstrapped_models/", MODEL_NAME, "_diff.Rds"))
}


