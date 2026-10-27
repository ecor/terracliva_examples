####
####  Monthly L-Moment / Regression Analysis and SPI on CHIRPS monthly precipitation
####  with readme (.readme.csv) files describing every output layer
####

library(terracliva)
library(magrittr)
library(raster)
library(lubridate)
library(data.table)
library(stringr)
library(terra)

wpath <- "/home/ecor/local/tasks/jrc/2026b/lm_and_spi"
wpath_data <- file.path(wpath, "data")
wpath_outcome <- file.path(wpath, "outcomes")

################
################ INPUT DATA
################

monthly_precipitation <- file.path(wpath_data, "chirps_monthly.tif") |> rast()
terra::time(monthly_precipitation) <- names(monthly_precipitation) %>%
  paste0("_01") %>% as.Date(format = "X%Y_%m_%d")

cores <- parallel::detectCores() - 4

################
################ README WRITERS
################

## Writes (and returns) a readme table describing the layers of a raster
## produced by lmapprast(). The CSV is saved next to the .tif file.
##
## lm          SpatRaster returned by lmapprast()
## x           SpatRaster of the analysed time series (used for the time range)
## var         variable short name (e.g. "tmax")
## unit        unit of the variable (e.g. "K")
## variable_description  long description of the variable
## filename    path of the .tif file; the CSV gets the same name with ".readme.csv"
## distrib     name of the fitted distribution (e.g. "pe3"), NULL if none
## product     data product name
## exc_label   tag used in the excess return-period layer names ({distrib}_{exc_label}_rt_%04d, {distrib}_{exc_label}_perc_rt_%04d)
## def_label   tag used in the deficit return-period layer names ({distrib}_{def_label}_rt_%04d, {distrib}_{def_label}_perc_rt_%04d)

write_lm_readme <- function(lm, x, var, unit, variable_description, filename, distrib = NULL,
                            product = "ERA5-Land", exc_label = "exc", def_label = "def") {
  
  l_moments_abs   <- sprintf("l_%d", 1:9)
  l_moment_ratios <- c("l_cv", "t", sprintf("t_%d", 1:9))
  regression_analysis_parameters <- c("pvalue", "coeff", "stdrerror", "rsquared", "senslope", "pvalue_mk")
  
  readme <- data.table(file = basename(filename), var = var, layer = names(lm))
  readme$product <- product
  readme$analysis_time_range_year_start <- min(year(terra::time(x)))
  readme$analysis_time_range_year_end   <- max(year(terra::time(x)))
  readme$variable_description <- variable_description
  readme$variable_unit <- unit
  readme$layer_unit <- "-"
  readme$month <- str_split(readme$layer, "_", n = 2) |> sapply(function(x) { x[[1]] })
  readme$index <- str_split(readme$layer, "_", n = 2) |> sapply(function(x) { x[[2]] })
  readme$is_l_moment <- readme$index %in% c(l_moments_abs, l_moment_ratios)
  is_l_moment_abs <- readme$index %in% c(l_moments_abs)
  readme$layer_unit[is_l_moment_abs] <- readme$variable_unit[is_l_moment_abs]
  readme$is_regression_analysis_parameter <- readme$index %in% regression_analysis_parameters
  ## return-period quantiles:
  ##   {distrib}_{def_label}_rt_%04d       deficit (variable unit)
  ##   {distrib}_{exc_label}_rt_%04d       excess  (variable unit)
  ##   {distrib}_{def_label}_perc_rt_%04d  deficit as percentage of the monthly mean (l_1)
  ##   {distrib}_{exc_label}_perc_rt_%04d  excess  as percentage of the monthly mean (l_1)
  readme$is_return_period_deficit <- FALSE
  readme$is_return_period_excess  <- FALSE
  readme$is_return_period_deficit_perc <- FALSE
  readme$is_return_period_excess_perc  <- FALSE
  readme$return_period_years <- NA_integer_
  if (!is.null(distrib)) {
    rp_patterns <- c(
      is_return_period_deficit      = sprintf("^%s_%s_rt_([0-9]+)$", distrib, def_label),
      is_return_period_excess       = sprintf("^%s_%s_rt_([0-9]+)$", distrib, exc_label),
      is_return_period_deficit_perc = sprintf("^%s_%s_perc_rt_([0-9]+)$", distrib, def_label),
      is_return_period_excess_perc  = sprintf("^%s_%s_perc_rt_([0-9]+)$", distrib, exc_label)
    )
    for (col in names(rp_patterns)) {
      is_rp <- str_detect(readme$index, rp_patterns[[col]])
      rp <- as.integer(str_match(readme$index, rp_patterns[[col]])[, 2])
      readme[[col]] <- is_rp
      readme$return_period_years[is_rp] <- rp[is_rp]
    }
  }
  is_return_period <- readme$is_return_period_deficit | readme$is_return_period_excess |
    readme$is_return_period_deficit_perc | readme$is_return_period_excess_perc
  ## any other layer is treated as a parameter of the fitted distribution (if any)
  readme$is_distribution_parameter <- !is.null(distrib) & !readme$is_l_moment & !readme$is_regression_analysis_parameter & !is_return_period
  readme$distribution <- "-"
  if (!is.null(distrib)) readme$distribution[readme$is_distribution_parameter | is_return_period] <- distrib
  readme$regression_analysis <- "-"
  readme$regression_analysis[readme$index %in% c("coeff", "stdrerror", "rsquared", "pvalue")] <- "Linear Model Analysis: variable vs time(years)"
  readme$regression_analysis[readme$index %in% c("senslope", "pvalue_mk")] <- "Mann-Kendall correlation test / Sen's slope estimation"
  readme$note <- "-"
  readme$note[readme$index %in% c("pvalue", "pvalue_mk")] <- "P-Value (if it is < 0.1) senslope or coeff are significant"
  readme$note[readme$index %in% c("senslope")] <- "Annual Sens-Slope Trend expressed in %s/year unit" |> sprintf(unit)
  readme$layer_unit[readme$index %in% c("senslope")] <- "%s/year" |> sprintf(readme$variable_unit[readme$index %in% c("senslope")])
  readme$note[readme$index %in% l_moments_abs] <- "Absolute L-Moment expressed in %s unit" |> sprintf(unit)
  readme$note[readme$index %in% l_moment_ratios] <- "L-Moment Ratio (dimensionless): l_cv=t ; t=l_2/l_1 ; t_i=l_i/l_2 if i>=3"
  if (!is.null(distrib)) {
    readme$note[readme$is_distribution_parameter] <- "Parameter of the '%s' distribution fitted with L-moments" |> sprintf(distrib)
    i_def <- readme$is_return_period_deficit
    i_exc <- readme$is_return_period_excess
    i_def_perc <- readme$is_return_period_deficit_perc
    i_exc_perc <- readme$is_return_period_excess_perc
    mm <- readme$month
    rp <- readme$return_period_years
    readme$note[i_def] <- "Deficit associated to a %d years return period : precipitation equal or lower than the mapped value occur with %04d year return period (absolute values)." |> sprintf(rp[i_def], rp[i_def])
    readme$note[i_exc] <- "Excess associated to a %d years return period : precipitation equal or greater than the mapped value occur with %04d year return period (absolute values)." |> sprintf(rp[i_exc], rp[i_exc])
    readme$note[i_def_perc] <- "Deficit associated to a %d years return period expressed as percentage of the monthly mean (l_1) : precipitation equal or lower than the mapped percentage of the monthly mean occur with %04d year return period (percenteage of %s_l_1)" |> sprintf(rp[i_def_perc], rp[i_def_perc], mm[i_def_perc])
    readme$note[i_exc_perc] <- "Excess associated to a %d years return period expressed as percentage of the monthly mean (l_1) : precipitation equal or greater than the mapped percentage of the monthly mean occur with %04d year return period(percenteage of %s_l_1)" |> sprintf(rp[i_exc_perc], rp[i_exc_perc], mm[i_exc_perc])
    
    readme$layer_unit[i_def] <- readme$variable_unit[i_def]
    readme$layer_unit[i_exc] <- readme$variable_unit[i_exc]
    readme$layer_unit[i_def_perc] <- "percentage"
    readme$layer_unit[i_exc_perc] <- "percentage"
  }
  
  filename_csv <- paste0(filename)
  extension(filename_csv) <- ".readme.csv"
  write.table(readme, file = filename_csv, row.names = FALSE, quote = FALSE, sep = ",")
  
  return(readme)
}

## Writes (and returns) a readme table describing the layers of a raster
## produced by spiapprast(). The CSV is saved next to the .tif file.
##
## Layer names produced by terracliva::spiapprast() / spicliva():
##   spi%02d_on_%04d_%02d   SPI value for each time step (year, month)
##   spi%02d_<param>        regression parameters (if summary_regress=TRUE):
##                          pvalue, coeff, stdrerror, rsquared, senslope, pvalue_mk
##
## spi         SpatRaster returned by spiapprast()
## x           SpatRaster of the analysed time series (used for the time range)
## var         variable short name (e.g. "precipitation")
## unit        unit of the input variable (e.g. "mm")
## variable_description  long description of the input variable
## filename    path of the .tif file; the CSV gets the same name with ".readme.csv"
## spi.scale   SPI aggregation scale in months (e.g. 3 for SPI-3)
## distrib     distribution fitted with L-moments (e.g. "pe3")
## timex_ref   dates used to fit the distribution (reference period)
## product     data product name
## pthres      tail probability used to clip SPI before the regression (spicliva default 1e-5)

write_spi_readme <- function(spi, x, var, unit, variable_description, filename,
                             spi.scale = attr(spi, "spi.scale"), distrib = "pe3",
                             timex_ref = terra::time(x), product = "CHIRPS-v3",
                             pthres = 10^-5) {
  
  if (is.null(spi.scale)) spi.scale <- 1
  regression_analysis_parameters <- c("pvalue", "coeff", "stdrerror", "rsquared", "senslope", "pvalue_mk")
  prefix <- sprintf("spi%02d", spi.scale)
  
  readme <- data.table(file = basename(filename), var = var, layer = names(spi))
  readme$product <- product
  readme$analysis_time_range_year_start <- min(year(terra::time(x)))
  readme$analysis_time_range_year_end   <- max(year(terra::time(x)))
  readme$reference_period_year_start <- min(year(timex_ref))
  readme$reference_period_year_end   <- max(year(timex_ref))
  readme$variable_description <- variable_description
  readme$variable_unit <- unit
  readme$spi_scale_months <- spi.scale
  readme$distribution <- distrib
  readme$layer_unit <- "-"
  
  ## time-step layers: spiSS_on_YYYY_MM
  ts_pattern <- sprintf("^%s_on_([0-9]{4})_([0-9]{2})$", prefix)
  readme$is_spi_value <- str_detect(readme$layer, ts_pattern)
  mt <- str_match(readme$layer, ts_pattern)
  readme$year  <- as.integer(mt[, 2])
  readme$month <- as.integer(mt[, 3])
  
  ## regression layers: spiSS_<param>
  readme$index <- str_replace(readme$layer, sprintf("^%s_", prefix), "")
  readme$index[readme$is_spi_value] <- "spi"
  readme$is_regression_analysis_parameter <- readme$index %in% regression_analysis_parameters
  
  readme$regression_analysis <- "-"
  readme$regression_analysis[readme$index %in% c("coeff", "stdrerror", "rsquared", "pvalue")] <- "Linear Model Analysis: SPI vs time(years)"
  readme$regression_analysis[readme$index %in% c("senslope", "pvalue_mk")] <- "Mann-Kendall correlation test / Sen's slope estimation"
  
  readme$note <- "-"
  i_spi <- readme$is_spi_value
  readme$note[i_spi] <- sprintf(
    "SPI-%d (dimensionless; standard normal) for %04d-%02d: %d-month cumulated %s ending in this month; transformed via '%s' distribution fitted with L-moments on %d-%d",
    spi.scale, readme$year[i_spi], readme$month[i_spi], spi.scale, tolower(variable_description),
    distrib, readme$reference_period_year_start[i_spi], readme$reference_period_year_end[i_spi])
  readme$note[readme$index %in% c("pvalue", "pvalue_mk")] <- "P-Value (if it is < 0.1) senslope or coeff are significant"
  readme$note[readme$index %in% "coeff"] <- "Linear trend of SPI expressed in SPI units/year"
  readme$note[readme$index %in% "senslope"] <- "Annual Sen's Slope trend of SPI expressed in SPI units/year"
  clip_note <- sprintf(" (SPI values clipped to +/-%.2f - tail probability %g - before regression)", -qnorm(pthres), pthres)
  i_reg <- readme$is_regression_analysis_parameter
  readme$note[i_reg] <- paste0(readme$note[i_reg], clip_note)
  readme$layer_unit[readme$index %in% c("coeff", "senslope")] <- "SPI/year"
  
  filename_csv <- paste0(filename)
  extension(filename_csv) <- ".readme.csv"
  write.table(readme, file = filename_csv, row.names = FALSE, quote = FALSE, sep = ",")
  
  return(readme)
}

################
################ LM ANALYSIS
################

var <- "precipitation"
distrib <- "pe3"
filename <- file.path(wpath_outcome, "lm_monthly_%s_variable.tif") |> sprintf(var)

if (file.exists(filename)) {
  lm_monthly_precipitation <- rast(filename)
} else {
  lm_monthly_precipitation <- lmapprast(monthly_precipitation, index = "monthly", distrib = distrib,
                                        summary_regress = TRUE, cores = cores,
                                        filename = filename, overwrite = TRUE)
}

lm_monthly_precipitation_readme <- write_lm_readme(lm_monthly_precipitation, monthly_precipitation,
                                                   var = var, unit = "mm",
                                                   variable_description = "Monthly Precipitation",
                                                   filename = filename, distrib = distrib,
                                                   product = "CHIRPS-v3")

################
################ SPI
################

spi.scale <- 3 ## put 1,3,6,12,24 alternatively
filename <- file.path(wpath_outcome, "spi_%02d_monthly_%s_variable.tif") |> sprintf(spi.scale, var)

## reference period for the distribution fitting
itimex_ref <- which(year(terra::time(monthly_precipitation)) %in% 1981:2010)
timex_ref <- terra::time(monthly_precipitation)[itimex_ref]
## NOTE: the following line overrides the 1981-2010 selection and uses the full record.
## Comment it out to use 1981-2010 as reference period.
timex_ref <- terra::time(monthly_precipitation)

if (file.exists(filename)) {
  spi_monthly_precipitation <- rast(filename)
} else {
  spi_monthly_precipitation <- spiapprast(x = monthly_precipitation, summary_regress = TRUE,
                                          spi.scale = spi.scale, filename = filename,
                                          distrib = distrib, cores = cores,
                                          timex_ref = timex_ref, overwrite = TRUE)
}

## spi.scale is passed explicitly: the attribute is lost when the raster is reloaded with rast()
spi_monthly_precipitation_readme <- write_spi_readme(spi_monthly_precipitation, monthly_precipitation,
                                                     var = var, unit = "mm",
                                                     variable_description = "Monthly Precipitation",
                                                     filename = filename, spi.scale = spi.scale,
                                                     distrib = distrib, timex_ref = timex_ref,
                                                     product = "CHIRPS-v3")


