library(raster);
library(tidyverse);
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

prismfiles <- list.files(prism_2pt5deg_path)
fileyears <- substr(prismfiles, 24, 27)
filemonths <- substr(prismfiles, 28, 29)

## read in files

prism_raster <- stack(paste0(prism_2pt5deg_path, "/", prismfiles[which(fileyears %in% prism_monthly_startyear:prism_monthly_endyear)]))

fileyears <- fileyears[which(fileyears %in% prism_monthly_startyear:prism_monthly_endyear)]
filemonths <- filemonths[which(fileyears %in% prism_monthly_startyear:prism_monthly_endyear)]


names(prism_raster) <- paste0(fileyears, "_", filemonths)

## save grid cell time series in data.frame (incl. coordinates for matching)
precip_df <- as.data.frame(prism_raster, xy = TRUE, na.rm = TRUE, centroids = TRUE,
                           long = TRUE) %>%
        dplyr::rename(date = layer, monthly_precip = value) %>%
        mutate(date = gsub("X", "",date)) %>%
        separate(date, into = c("year", "month"), sep = "_") %>%
        mutate(centroid_coord = paste0(x, "_", y),
               grid_res = xres(prism_raster))

saveRDS(precip_df, file = monthly_prism_2pt5_file)


