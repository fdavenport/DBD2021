library(raster); library(velox);
library(dplyr); library(tidyr);
source("./file_paths.R")
source("./func.R")
source("./parameters.R")

## -----------------------------------------------------------------------------
## READ IN STATE BOUNDARIES
print("Reading in boundary files...")
statebounds <- shapefile(state_boundary_file)

## MISSOURI BOUNDARY
MO_bounds <- subset(statebounds, NAME == "Missouri")
## MISSISSIPPI RIVER BOUNDARY
MIR_bounds <- shapefile(um_basin_file)
MIR_noMO_bounds <- MIR_bounds - MO_bounds

## MISSOURI RIVER BOUNDARY
MOR_bounds <- shapefile(mo_basin_file)
MOR_noMO_bounds <- MOR_bounds - MO_bounds

MOR_bounds$Name <- "Missouri_River"
MIR_bounds$Name <- "Mississippi_River"
MIR_noMO_bounds$Name <- "Mississippi_noMO"
MOR_noMO_bounds$Name <- "Missouri_noMO"
basin_bounds <- do.call(rbind,
                            list(MOR_bounds, MIR_bounds, MIR_noMO_bounds, MOR_noMO_bounds))
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
## SEQUENCE USES FIRST DAY OF MONTH BUT DATA IS FOR CUMULATIVE MONTHLY AMOUNT
month_dates <- seq.Date(as.Date(paste0(prism_monthly_startyear, "-01-01")),
                        as.Date(paste0(prism_monthly_endyear, "-12-01")), by = "month")

## -----------------------------------------------------------------------------
## CALCULATE STATE AVERAGES FOR EACH MONTH
prismVelox <- velox(prismRaster)
print("Extracting monthly PPT for each area...")
extractedPPT <- prismVelox$extract(sp = basin_bounds, small = TRUE)
## calculate watershed mean (done outside of extract to allow for NA remove)
monthlyPPT_list <- lapply(extractedPPT, function(x) colMeans(x, na.rm = TRUE))

## assign dates
monthlyPrecipDF <- as.data.frame(t(matrix(unlist(monthlyPPT_list),
                                          nrow = length(month_dates)))) %>%
    setNames(., as.character(strftime(month_dates, "%Y-%m"))) %>%
    mutate(NAME = basin_bounds$Name) %>%
    gather(date, monthly_precip, as.character(strftime(month_dates, "%Y-%m"))) %>%
    separate(date, into = c("year", "month"), sep = "-") %>%
    mutate(year = as.numeric(year))

## save data
saveRDS(monthlyPrecipDF, watershed_precip_file)
