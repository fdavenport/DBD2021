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
## READ MONTHLY PRISM DATA
prismfiles <- list.files(prism_monthly_path, pattern = "\\.bil") ## use \\ to include "." in pattern,
prismfiles <- prismfiles[seq(1, length(prismfiles)-1, 2)] ## remove .aux.xml files that also get included
## subset for files starting after certain date
fileyears <- substr(prismfiles, 24, 27)
prismfiles <- prismfiles[which(fileyears %in% prism_monthly_startyear:prism_monthly_endyear)]
prismfiles <- paste0(prism_monthly_path, "/", prismfiles)
print("Reading in PRISM monthly files...")
prismRaster <- stack(prismfiles)

## SEQUENCE OF DATES FOR ALL MONTHS
## NOTE THAT THIS USES THE FIRST DAY OF MONTH BUT DATA IS FOR CUMULATIVE MONTHLY AMOUNT
month_dates <- seq.Date(as.Date(paste0(prism_monthly_startyear, "-01-01")),
                        as.Date(paste0(prism_monthly_endyear, "-12-01")), by = "month")

## -----------------------------------------------------------------------------
## CALCULATE STATE AVERAGES FOR EACH MONTH
prismVelox <- velox(prismRaster)
print("Extracting monthly PPT for each state...")
extractedPPT <- prismVelox$extract(sp = statebounds, small = TRUE)
## check for states not covered by PRISM data and remove
states_missing_ppt <- which(sapply(extractedPPT, is.null))
extractedPPT <- extractedPPT[-states_missingPPT]
## calculate state mean (done outside of extract to allow for NA remove)
monthlyPPT_list <- lapply(extractedPPT, function(x) colMeans(x, na.rm = TRUE))

## assign dates
monthlyPrecipDF <- as.data.frame(t(matrix(unlist(monthlyPPT_list),
                                          nrow = length(month_dates)))) %>%
    setNames(., as.character(strftime(month_dates, "%Y-%m"))) %>%
    mutate(STATESFP = statebounds$STATEFP[-states_missing_ppt]) %>%
    mutate(STATENAME = statebounds$NAME[-states_missing_ppt]) %>%
    gather(date, monthly_precip, as.character(strftime(month_dates, "%Y-%m"))) %>%
    separate(date, into = c("year", "month"), sep = "-") %>%
    mutate(year = as.numeric(year))

## save data
saveRDS(monthlyPrecipDF, monthly_precip_file)
