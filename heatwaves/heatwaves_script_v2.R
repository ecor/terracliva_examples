## ---------------------------------------------------------------------------
## Install missing packages
## ---------------------------------------------------------------------------
cran_pkgs <- c("terra", "raster", "data.table", "stringr", "leaflet")  ## leaflet is needed by terra::plet()
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

if (!requireNamespace("terracliva", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("ecor/terracliva")
}

## ---------------------------------------------------------------------------
## Load packages
## ---------------------------------------------------------------------------
library(terracliva)
library(terra)
library(parallel)   ## base R, always available
library(raster)
library(data.table)
library(stringr)

## ---------------------------------------------------------------------------
## Paths and settings
## ---------------------------------------------------------------------------
wpath <- "/home/ecor/local/tasks/jrc/2026b/terracliva_examples/heatwaves"
wpath_data <- file.path(wpath, "data")
wpath_outcome <- file.path(wpath, "outcomes")
cores <- parallel::detectCores() - 4

## ---------------------------------------------------------------------------
## Input data
## ---------------------------------------------------------------------------
tmax_dataset_daily <- wpath_data |> list.files(full.names = TRUE, pattern = "tmax.nc") |> rast() 
terra::time(tmax_dataset_daily) <- terra::time(tmax_dataset_daily) |> as.Date()

var <- "tmax"

filename <- file.path(wpath_outcome, "hwmid_%s_variable.tif") |> sprintf(var)

## ---------------------------------------------------------------------------
## HWMId computation (skipped if the output already exists)
## ---------------------------------------------------------------------------
if (file.exists(filename)) {
  o_hw_regress <- rast(filename)
} else {
  o_hw_regress <- hwmidapprast(tmax_dataset_daily, summary_regress = TRUE, cores = cores,
                               filename = filename, overwrite = TRUE)
}

## ---------------------------------------------------------------------------
## write_hwmid_readme
## Writes a <filename>.readme.csv file describing, layer by layer, the output
## of hwmidapprast() / cwmidapprast() (with summary_regress=TRUE or FALSE).
##
## Layers produced by hwmidcliva():
##   - "YYYY" (one per year)    : annual value of the HWMId/CWMId index
##   - "<pfx>_pvalue", "<pfx>_coeff", "<pfx>_stdrerror", "<pfx>_rsquared",
##     "<pfx>_senslope", "<pfx>_pvalue_mk"  : from terracliva::regress()
##   - "<pfx>_n years with index equal or greater than <thres>"
##   where <pfx> = "hwmid" (heat wave) or "cwmid" (cold wave).
##
## Arguments:
##   hw            SpatRaster returned by hwmidapprast()/cwmidapprast()
##   x             daily input SpatRaster (e.g. tmax)
##   var, unit, variable_description  description of the input variable
##   filename      raster file of the output (the readme is written next to it)
##   start_month   first month of the analysis "year"; "default" mirrors the
##                 hwmidapprast() rule (1 northern hemisphere, 7 southern)
##   timex_ref     dates of the reference period (default: the whole of x,
##                 as in hwmidcliva())
##   hwmid_thres   threshold used to count the years (default 4)
##   signif        significance level used in regress() (default 0.1)
##   cold          TRUE if the output comes from cwmidapprast()
##   product       name of the source dataset
## ---------------------------------------------------------------------------
write_hwmid_readme <- function(hw, x, var, unit, variable_description, filename,
                               start_month = "default",
                               timex_ref = terra::time(x),
                               hwmid_thres = 4, signif = 0.1, cold = FALSE,
                               product = "-") {
  
  ## --- first month of the analysis year: same rule as hwmidapprast() ------
  if (start_month == "default") {
    ee <- x[[1]] |> terra::project("epsg:4326") |> terra::ext()
    start_month <- if ((ee$ymin + ee$ymax) / 2 >= 0) 1 else 7
  }
  start_month <- as.integer(start_month)
  
  ## --- index prefix and description --------------------------------------
  prefix     <- if (cold) "cwmid" else "hwmid"
  index_name <- if (cold) "CWMId (Cold Wave Magnitude Index daily)" else
    "HWMId (Heat Wave Magnitude Index daily)"
  regression_analysis_parameters <- c("pvalue", "coeff", "stdrerror", "rsquared",
                                      "senslope", "pvalue_mk")
  count_layer <- sprintf("n years with index equal or greater than %d", hwmid_thres)
  
  ## --- base table ---------------------------------------------------------
  readme <- data.table(file = basename(filename), var = var, layer = names(hw))
  readme$product <- product
  readme$analysis_time_range_year_start <- min(year(terra::time(x)))
  readme$analysis_time_range_year_end   <- max(year(terra::time(x)))
  readme$reference_period_year_start    <- min(year(timex_ref))
  readme$reference_period_year_end      <- max(year(timex_ref))
  readme$variable_description <- variable_description
  readme$variable_unit        <- unit
  readme$index_name           <- index_name
  readme$start_month          <- start_month
  readme$layer_unit           <- "-"
  
  ## --- annual layers: "YYYY" (optional leading "X" if names were mangled) -
  ts_pattern <- "^X?([0-9]{4})$"
  readme$is_index_value <- str_detect(readme$layer, ts_pattern)
  readme$year <- as.integer(str_match(readme$layer, ts_pattern)[, 2])
  
  ## time span covered by each analysis "year" (with start_month != 1 it
  ## spans two calendar years, e.g. July YYYY - June YYYY+1)
  readme$period_start <- NA_character_
  readme$period_end   <- NA_character_
  i_idx <- readme$is_index_value
  ps <- as.Date(sprintf("%04d-%02d-01", readme$year[i_idx], start_month))
  readme$period_start[i_idx] <- format(ps)
  readme$period_end[i_idx]   <- format(as.Date(sapply(ps, function(d)
    format(seq(d, by = "1 year", length.out = 2)[2] - 1))))
  
  ## --- summary layers: "<prefix>_<param>" ---------------------------------
  readme$index <- str_replace(readme$layer, sprintf("^%s_", prefix), "")
  readme$index[i_idx] <- prefix
  readme$is_regression_analysis_parameter <- readme$index %in% regression_analysis_parameters
  readme$is_count <- readme$index == count_layer
  
  readme$regression_analysis <- "-"
  readme$regression_analysis[readme$index %in% c("coeff", "stdrerror", "rsquared", "pvalue")] <-
    sprintf("Linear Model Analysis: %s vs time(years)", toupper(prefix))
  readme$regression_analysis[readme$index %in% c("senslope", "pvalue_mk")] <-
    "Mann-Kendall correlation test / Sen's slope estimation"
  
  ## --- notes and units ----------------------------------------------------
  readme$note <- "-"
  readme$note[i_idx] <- sprintf(
    "%s for the year starting in %s (%s to %s), computed on daily %s; thresholds from reference period %d-%d",
    index_name, month.name[start_month], readme$period_start[i_idx], readme$period_end[i_idx],
    tolower(variable_description),
    readme$reference_period_year_start[i_idx], readme$reference_period_year_end[i_idx])
  readme$note[readme$index %in% c("pvalue", "pvalue_mk")] <-
    sprintf("P-Value (if it is < %g) senslope or coeff are significant", signif)
  readme$note[readme$index %in% "coeff"] <-
    sprintf("Linear trend of %s expressed in index units/year", toupper(prefix))
  readme$note[readme$index %in% "senslope"] <-
    sprintf("Annual Sen's Slope trend of %s expressed in index units/year", toupper(prefix))
  readme$note[readme$index %in% "stdrerror"] <- "Standard error of the linear trend coefficient"
  readme$note[readme$index %in% "rsquared"]  <- "R squared of the linear model"
  
  ## in regress(), coeff, stdrerror, rsquared and senslope are set to NA where
  ## not significant
  i_masked <- readme$index %in% c("coeff", "stdrerror", "rsquared", "senslope")
  readme$note[i_masked] <- paste0(readme$note[i_masked],
                                  sprintf(" (set to NA where the corresponding p-value > %g)", signif))
  
  readme$note[readme$is_count] <- sprintf(
    "Number of years in %d-%d with %s >= %d",
    readme$analysis_time_range_year_start[readme$is_count],
    readme$analysis_time_range_year_end[readme$is_count],
    toupper(prefix), hwmid_thres)
  
  readme$layer_unit[readme$index %in% c("coeff", "senslope", "stdrerror")] <-
    sprintf("%s/year", toupper(prefix))
  readme$layer_unit[readme$is_count] <- "years"
  
  ## --- write --------------------------------------------------------------
  filename_csv <- paste0(filename)
  raster::extension(filename_csv) <- ".readme.csv"
  write.table(readme, file = filename_csv, row.names = FALSE, quote = FALSE, sep = ",")
  
  return(readme)
}

## ---------------------------------------------------------------------------
## Write readme
## start_month must be given explicitly if hwmidapprast() was not called
## with start_month = "default"
## ---------------------------------------------------------------------------
o_hw_regress_readme <- write_hwmid_readme(o_hw_regress, tmax_dataset_daily,
                                          var = var, unit = "degC",
                                          variable_description = "Daily Maximum Temperature",
                                          filename = filename,
                                          product = "ERA5Land")

## ---------------------------------------------------------------------------
## Interactive map
## ---------------------------------------------------------------------------
plet(o_hw_regress$hwmid_senslope, tiles = "OpenTopoMap")