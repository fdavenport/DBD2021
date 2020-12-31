library(raster)
library(velox)
library(dplyr)
library(tidyr)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")


## -----------------------------------------------------------------------------
## READ IN STATE BOUNDARIES
print("Reading in state boundary file...")
statebounds <- shapefile(state_boundary_file)

## -----------------------------------------------------------------------------
## READ IN DAILY PRISM DATA
print("Reading in PRISM files...")

daily_ppt_list <- vector('list', length = length(prism_daily_years))
for(i in seq_along(prism_daily_years)){
    print(prism_daily_years[i])
    prismfiles <- list.files(paste0(prism_daily_path, "/", prism_daily_years[i]),
                             pattern = "\\.bil$", full.names = TRUE)

    prism_raster <- stack(prismfiles)
    prism_velox <- velox(prism_raster)

    n_days <- dim(prism_raster)[3]
    yeardates <- seq.Date(as.Date(paste0(prism_daily_years[i], "-01-01")), by = "day",
                          length.out = n_days)

    ## CALCULATE DAILY AVERAGE FOR EACH COUNTY
    ##print("Extracting daily PPT for each county...")
    extracted_ppt <- prism_velox$extract(sp = statebounds)
    ## check for states not covered by PRISM data and remove
    states_missing_ppt<- which(sapply(extracted_ppt, is.null))
    extracted_ppt <- extracted_ppt[-states_missing_ppt]
    ## calculate mean value for each state
    extracted_ppt <- lapply(extracted_ppt, function(x) colMeans(x, na.rm = TRUE))

    daily_ppt_df <- as.data.frame(t(matrix(unlist(extracted_ppt),
                                          nrow = length(yeardates)))) %>%
    setNames(., as.character(strftime(yeardates, "%Y-%m-%d"))) %>%
    mutate(STATESFP = statebounds$STATEFP[-states_missing_ppt]) %>%
    mutate(STATENAME = statebounds$NAME[-states_missing_ppt]) %>%
    gather(date, daily_precip, as.character(strftime(yeardates, "%Y-%m-%d")))

    daily_ppt_list[[i]] <- daily_ppt_df
}

daily_ppt <- bind_rows(daily_ppt_list)

saveRDS(daily_ppt, daily_precip_file)
