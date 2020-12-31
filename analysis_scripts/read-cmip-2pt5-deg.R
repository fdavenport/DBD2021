library(raster)
library(ncdf4)
library(dplyr)
library(tidyr)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

## READ IN COMMAND LINE ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
ARRAY_ID <- as.numeric(args[1]) ## job number within array
print(paste("ARRAY_ID = ", ARRAY_ID))
FORCING <- args[2]

if(FORCING == "historical"){
    filepath <- cmip_hist_path
} else if (FORCING == "rcp85"){
    filepath <- cmip_rcp85_path
} else if (FORCING == "rcp26"){
    filepath <- cmip_rcp26_path
} else {
    stop("Incorrect forcing specified")
}

MODEL_NAME <- list.files(filepath)[ARRAY_ID]
RUN_NAMES <- list.files(paste0(filepath, "/", MODEL_NAME))

## -----------------------------------------------------------------------------
## READ IN GRIDDED TIME SERIES FOR EACH RUN

for (i in 1:length(RUN_NAMES)){
    print(paste("Run:", RUN_NAMES[i]))

    currentfiles <- list.files(paste0(filepath, "/", MODEL_NAME,
                                      "/", RUN_NAMES[i]))

    print("Files:")
    print(currentfiles)
    ## CHECK STACK FOR MODEL WITH MULTIPLE FILES
    currentdata <- stack(paste0(filepath, "/", MODEL_NAME, "/",
                                RUN_NAMES[i], "/", currentfiles))
    print(paste("crs:", crs(currentdata)))
    print(paste("x-resolution:", xres(currentdata)))
    print(paste("y-resolution:", yres(currentdata)))
    print(paste("Number of Dates:", dim(currentdata)[3]))
    print("Date Format:")
    print(head(names(currentdata)))

    ## -------------------------------------------------------------------------
    ## CROP RASTER TO PRISM EXTENT (read in example file)
    prism_mask <- raster(list.files(path = prism_2pt5deg_path,
                                    full.names = TRUE)[1])

    currentdata <- mask(currentdata, prism_mask)

     ## create date sequence from file names
    start_end_dates <- lapply(currentfiles,
                              function(x) strsplit(strsplit(x, "_")[[1]][6], "-")[[1]])
    date_sequence <- do.call("c",
                             lapply(start_end_dates,
                                    function(x) seq.Date(as.Date(paste0(x[1], "15"), "%Y%m%d"),
                                                         as.Date(paste0(x[2], "15"), "%Y%m%d"),
                                                         by = "month")))
    ## CHECK FOR DUPLICATED DATES AND REMOVE FROM DATA!
    dup_date <- which(duplicated(date_sequence))
    if(length(dup_date) > 0){
        date_sequence <- date_sequence[-dup_date]
        currentdata <- dropLayer(currentdata, dup_date)
        print("OVERLAPPING DATES DETECTED AND REMOVED!!")
    }
    ## CHECK THAT DATE SEQUENCE AND NUMBER OF RASTER LAYERS MATCH
    if(nlayers(currentdata) != length(date_sequence)) {
        warning("Number of raster layers and number of dates do not match")
    }
    try(names(currentdata) <- date_sequence)

    precip_df <- as.data.frame(currentdata, xy = TRUE, na.rm = TRUE, centroids = TRUE,
                               long = TRUE) %>%
        dplyr::rename(date = layer, monthly_precip = value) %>%
        mutate(date = as.Date(gsub("X", "", gsub("\\.", "-", date))),
               year = strftime(date, "%Y"),
               month = strftime(date, "%m"),
               centroid_coord = paste0(x, "_", y),
               grid_res = xres(currentdata),
               FORCING = FORCING,
               MODEL = MODEL_NAME,
               RUN = RUN_NAMES[i]) %>%
        dplyr::select(-date)

    saveRDS(precip_df, file = paste0(monthly_cmip_path, MODEL_NAME, "_",
                                     FORCING, "_", RUN_NAMES[i], "_2pt5deg_precip.Rds"))
}
