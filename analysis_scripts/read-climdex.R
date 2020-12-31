library(tidyverse)
library(ncdf4)
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

rx5day <- nc_open(climdex_obs_rx5day_file)

lons <- ncvar_get(rx5day, "lon")
lats <- ncvar_get(rx5day, "lat")
time <- as.Date(ncvar_get(rx5day, "time"), origin = "1850-01-15")

monthly_rx5day <- expand.grid(lons, lats, time) %>%
    dplyr::rename(lon = Var1, lat = Var2, time = Var3) %>%
    mutate(rx5day = ncvar_get(rx5day, "Rx5day"),
           year = strftime(time, "%Y"),
           month = strftime(time, "%m"))

nc_close(rx5day)

monthly_rx5day <- monthly_rx5day %>%
    subset(lon >= 235 & lon <= 292.5 & lat >= 26 & lat <= 49) %>%
    na.omit()

saveRDS(monthly_rx5day, rx5day_obs_file)

## -----------------------------------------------------------------------------

forcings <- c("historical", "rcp26", "rcp85")

for(i in forcings){
    files <- list.files(paste0(climdex_cmip_path, i), pattern="rx5day",
                        full.names = TRUE)

    for(j in seq_along(files)){
        print(j)
        current_file <- nc_open(files[j])
        model <- strsplit(files[j], "_")[[1]][4]
        runname <- strsplit(files[j], "_")[[1]][6]

        ## only read lon/lat over CONUS (approximate) to improve speed
        lons <- ncvar_get(current_file, "lon")
        ## sometimes there are two closest values, so take the first
        lon_start <- which(abs(lons - 232) == min(abs(lons-232)))[1]
        lon_end <- which(abs(lons-294) == min(abs(lons-294)))[1]
        lons <- lons[lon_start:lon_end]

        lats <- ncvar_get(current_file, "lat")
        lat_start <- which(abs(lats - 24) == min(abs(lats-24)))[1]
        lat_end <- which(abs(lats-50) == min(abs(lats-50)))[1]
        lats <- lats[lat_start:lat_end]

        start_end_dates <- strsplit(strsplit(files[j], "_")[[1]][7], "-")[[1]]
        date_sequence <-  seq.Date(as.Date(paste0(start_end_dates[1], "15"),
                                           "%Y%m%d"),
                                   as.Date(paste0(start_end_dates[2], "15"),
                                           "%Y%m%d"),
                                   by = "month")

        current_data <- expand.grid(lons, lats, date_sequence) %>%
            dplyr::rename(lon = Var1, lat = Var2, time = Var3) %>%
            mutate(rx5day = ncvar_get(current_file, "rx5day",
                                      start=c(lon_start, lat_start, 1),
                                      count=c(length(lons), length(lats), -1)),
                   MODEL = model,
                   RUN = runname,
                   FORCING = i,
                   year = strftime(time, "%Y"),
                   month = strftime(time, "%m"))

        saveRDS(current_data, paste0(rx5day_cmip_path,
                                     model, "_", i, "_", runname,
                                     "_2pt5deg_rx5day.Rds"))
    }
}


