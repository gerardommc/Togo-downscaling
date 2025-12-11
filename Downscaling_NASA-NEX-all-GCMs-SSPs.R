# ==========================================
# NEX-GDDP-CMIP6 Downscaling + Bias Correction
# using pre-calculated DOY biases (baseline 1979–2014)
# and a custom 10 km grid
# For one GCM at a time (example: GFDL-ESM4)
# This script was developed by Sandra Milena Castaño Quintero, https://github.com/Salamanadra
# ==========================================

library(terra)
library(lubridate)
library(dplyr)
library(fields) # for Tps (Thin Plate Spline)

##### DIRECTORY PATHS (ADJUST THESE) ---------------------------------------

# Clipped TIF files from NEX-GDDP-CMIP6 (Historical + Future) within AOI
# Expected file structure: VAR_MODEL_SCENARIO_DATE.tif
dir_nex <- "PATH_TO_NEX_GCM_TIFFS_INPUT/"

# NetCDF files with pre-calculated DOY bias (coarse resolution)
# Example expected filename: baseline_hurs_GFDL-ESM4_1979_2014_Bias_DOY.nc
dir_bias <- "PATH_TO_DOY_BIAS_NETCDF_INPUT/"

# Custom 10 km grid (reference template)
grid10km_path <- "PATH_TO_CUSTOM_GRID_TEMPLATE/grid10km.nc"

# Output directory for corrected NetCDF files
dir_out <- "PATH_TO_DOWNSCALED_NETCDF_OUTPUT/"

# GCM to process (change this for each model run)
model <- "GFDL-ESM4" # change to "IPSL-CM6A-LR", etc.

# Baseline period used for bias calculation
baseline_start <- as.Date("1979-01-01")
baseline_end   <- as.Date("2014-12-31")

# Epsilon value for logarithmic correction (to avoid log(0))
eps_pr   <- 0.1  # mm/day (Precipitation)
eps_hurs <- 0.01 # Humidity (adjust if input is 0–1 instead of 0–100)

##### LOAD 10 KM GRID TEMPLATE -------------------------------------

grid10km <- rast(grid10km_path)
crs(grid10km) <- "EPSG:4326" # Ensure correct CRS

# Assign placeholder values if the template raster has no values (to avoid warnings)
if (!hasValues(grid10km)) grid10km[] <- 1L

##### AUXILIARY FUNCTIONS ----------------------------------

# Day Of Year (DOY) calculation
doy_vec <- function(dates) yday(dates)

# Thin Plate Spline (TPS) interpolation from coarse bias layer to fine grid
tps_to_fine <- function(bias_layer, template_fine) {
  # Extract valid values from the coarse bias layer
  vals <- values(bias_layer, mat = FALSE)
  ok <- !is.na(vals)
  
  if (!any(ok)) {
    warning("All cells in bias_layer are NA; returning NA in the fine layer.")
    r_out <- template_fine
    values(r_out) <- NA
    return(r_out)
  }
  
  # Coordinates of valid cells
  cells_ok <- which(ok)
  pts_ok   <- terra::xyFromCell(bias_layer, cells_ok)
  vals_ok  <- vals[ok]
  
  # Fit Thin Plate Spline model
  tps_fit <- fields::Tps(pts_ok, vals_ok)
  
  # Coordinates of all cells in the fine grid
  fine_cells <- seq_len(ncell(template_fine))
  fine_xy    <- terra::xyFromCell(template_fine, fine_cells)
  
  # Predict values onto the fine grid
  pred_vals <- predict(tps_fit, fine_xy)
  
  r_out <- template_fine
  values(r_out) <- pred_vals
  r_out
}

# Apply TPS to every layer in a bias stack
tps_stack_to_fine <- function(bias_stack, template_fine) {
  nb <- nlyr(bias_stack)
  out_list <- vector("list", nb)
  for (i in 1:nb) {
    cat("TPS layer", i, "of", nb, "\n")
    out_list[[i]] <- tps_to_fine(bias_stack[[i]], template_fine)
  }
  rast(out_list)
}

# Fine grid reference: ALWAYS the 10 km grid template
load_template_fine <- function(var) {
  grid10km
}

# ---------------------------------
# Load NEX baseline (historical) for a given variable/model
# ---------------------------------
load_nex_baseline <- function(var, model) {
  files_nex <- list.files(
    dir_nex,
    pattern = sprintf("^%s_%s_historical_.*\\.tif$", var, model),
    full.names = TRUE
  )
  if (length(files_nex) == 0) {
    stop("No NEX baseline files found for: ", var, " - ", model)
  }
  
  r_all <- rast(files_nex)
  
  # NASA UNITS CONVERSION
  if (var == "pr") {
    # pr: kg m^-2 s^-1 → mm/day
    r_all <- r_all * 86400
  } else if (var %in% c("tas", "tasmin", "tasmax")) {
    # temps: Kelvin → °C
    r_all <- r_all - 273.15
  }
  # hurs remains in %, no change
  
  # Dates 1950–2014 (adjust if your series starts/ends differently)
  dates_all <- seq(as.Date("1950-01-01"), as.Date("2014-12-31"), by = "day")
  dates_all <- dates_all[seq_len(nlyr(r_all))]
  time(r_all) <- dates_all
  
  # Select only the baseline period (1979-2014)
  sel <- time(r_all) >= baseline_start & time(r_all) <= baseline_end
  r_all[[sel]]
}

# ---------------------------------
# Load DOY bias (coarse) and interpolate to 10 km grid
# ---------------------------------
load_bias_doy_fine <- function(var, model, template_fine) {
  # Adjust this pattern to match your actual Bias NetCDF filenames
  bias_file <- file.path(
    dir_bias,
    sprintf("baseline_%s_%s_1979_2014_Bias_DOY.nc", var, model)
  )
  
  if (!file.exists(bias_file)) {
    stop("Bias file not found: ", bias_file)
  }
  
  bias_coarse <- rast(bias_file)
  cat("DOY Bias loaded from: ", bias_file, "\n")
  
  # TPS interpolation to the 10 km grid
  bias_fine <- tps_stack_to_fine(bias_coarse, template_fine)
  bias_fine
}

# ---------------------------------
# Load NEX future (ssp245, ssp585) and resample to 10 km grid
# ---------------------------------
load_nex_future_fine <- function(var, model, scenario, template_fine) {
  pat <- sprintf("^%s_%s_%s_.*\\.tif$", var, model, scenario)
  files_fut <- list.files(dir_nex, pattern = pat, full.names = TRUE)
  
  if (length(files_fut) == 0) {
    stop("No NEX future files found for: ", var, " - ", model, " - ", scenario)
  }
  
  r_fut <- rast(files_fut)
  
  # NASA UNITS CONVERSION
  if (var == "pr") {
    # pr: kg m^-2 s^-1 → mm/day
    r_fut <- r_fut * 86400
  } else if (var %in% c("tas", "tasmin", "tasmax")) {
    # temps: Kelvin → °C
    r_fut <- r_fut - 273.15
  }
  # hurs in %, no change
  
  # Future dates (adjust if your downloaded range is smaller)
  dates_fut <- seq(as.Date("2015-01-01"), as.Date("2100-12-31"), by = "day")
  dates_fut <- dates_fut[seq_len(nlyr(r_fut))]
  time(r_fut) <- dates_fut
  
  # Resample to the fine grid (your grid10km)
  resample(r_fut, template_fine, method = "bilinear")
}

##### APPLY ADDITIVE BIAS CORRECTION (tas, tasmin, tasmax) --------
# Corrected = NEX_fine - Bias_DOY_fine

apply_bias_add_baseline <- function(var, model) {
  # 1) Load NEX baseline (coarse)
  r_nex <- load_nex_baseline(var, model)
  dates_b <- time(r_nex)
  
  # 2) Fine grid (10 km)
  template_fine <- load_template_fine(var)
  
  # 3) Resample NEX to fine grid
  r_nex_fine <- resample(r_nex, template_fine, method = "bilinear")
  time(r_nex_fine) <- dates_b # Re-assign dates
  
  # 4) Info for the loop
  doy_b <- yday(dates_b)
  nb_b  <- nlyr(r_nex_fine)
  
  # 5) DOY Bias (fine)
  bias_fine <- load_bias_doy_fine(var, model, template_fine)
  
  out_list <- vector("list", nb_b)
  for (i in 1:nb_b) {
    d <- doy_b[i]
    bias_layer <- bias_fine[[d]]
    nex_layer  <- r_nex_fine[[i]]
    # Additive correction: Corrected = NEX_fine - Bias
    out_list[[i]] <- nex_layer - bias_layer
  }
  
  r_corr <- rast(out_list)
  time(r_corr) <- dates_b
  r_corr
}


apply_bias_add_future <- function(var, model, scenario) {
  template_fine <- load_template_fine(var)
  bias_fine     <- load_bias_doy_fine(var, model, template_fine)
  
  r_fut_fine <- load_nex_future_fine(var, model, scenario, template_fine)
  
  dates_f <- time(r_fut_fine)
  doy_f   <- yday(dates_f)
  nb_f    <- nlyr(r_fut_fine)
  
  out_list <- vector("list", nb_f)
  for (i in 1:nb_f) {
    d <- doy_f[i]
    bias_layer <- bias_fine[[d]]
    nex_layer  <- r_fut_fine[[i]]
    out_list[[i]] <- nex_layer - bias_layer
  }
  
  r_corr <- rast(out_list)
  time(r_corr) <- dates_f
  r_corr
}

##### APPLY LOGARITHMIC BIAS CORRECTION (pr, hurs) ----------------------
# Corrected = exp( log(NEX + eps) - log(Bias) ) - eps

apply_bias_log_baseline <- function(var, model) {
  eps <- if (var == "pr") eps_pr else if (var == "hurs") eps_hurs else 0
  
  r_nex <- load_nex_baseline(var, model)
  dates_b <- time(r_nex)
  
  template_fine <- load_template_fine(var)
  r_nex_fine <- resample(r_nex, template_fine, method = "bilinear")
  time(r_nex_fine) <- dates_b
  
  doy_b <- yday(dates_b)
  nb_b  <- nlyr(r_nex_fine)
  
  # Bias for log correction is defined as log(NEX/OBS) = log(NEX) - log(OBS)
  # So, the final correction is: log(OBS) = log(NEX) - log(Bias)
  log_bias_fine <- load_bias_doy_fine(var, model, template_fine)
  
  out_list <- vector("list", nb_b)
  for (i in 1:nb_b) {
    d <- doy_b[i]
    nex_layer <- r_nex_fine[[i]]
    
    # Clip negative values before log transform
    vals <- values(nex_layer)
    vals[vals < 0] <- 0
    values(nex_layer) <- vals
    
    log_nex <- log(nex_layer + eps)
    log_bias_day <- log_bias_fine[[d]]
    
    log_corr <- log_nex - log_bias_day
    var_corr <- exp(log_corr) - eps # Inverse transformation
    
    # Final clip of negative results
    vals_corr <- values(var_corr)
    vals_corr[vals_corr < 0] <- 0
    values(var_corr) <- vals_corr
    
    out_list[[i]] <- var_corr
  }
  
  r_corr <- rast(out_list)
  time(r_corr) <- dates_b
  r_corr
}


apply_bias_log_future <- function(var, model, scenario) {
  eps <- if (var == "pr") eps_pr else if (var == "hurs") eps_hurs else 0
  
  template_fine <- load_template_fine(var)
  log_bias_fine <- load_bias_doy_fine(var, model, template_fine)
  
  r_fut_fine <- load_nex_future_fine(var, model, scenario, template_fine)
  
  dates_f <- time(r_fut_fine)
  doy_f   <- yday(dates_f)
  nb_f    <- nlyr(r_fut_fine)
  
  out_list <- vector("list", nb_f)
  for (i in 1:nb_f) {
    d <- doy_f[i]
    nex_layer <- r_fut_fine[[i]]
    
    # Clip negative values before log transform
    vals <- values(nex_layer)
    vals[vals < 0] <- 0
    values(nex_layer) <- vals
    
    log_nex <- log(nex_layer + eps)
    log_bias_day <- log_bias_fine[[d]]
    
    log_corr <- log_nex - log_bias_day
    var_corr <- exp(log_corr) - eps
    
    # Final clip of negative results
    vals_corr <- values(var_corr)
    vals_corr[vals_corr < 0] <- 0
    values(var_corr) <- vals_corr
    
    out_list[[i]] <- var_corr
  }
  
  r_corr <- rast(out_list)
  time(r_corr) <- dates_f
  r_corr
}


##### APPLICATION EXAMPLES FOR GFDL-ESM4 (Uncomment to Run) -------------------

# 1) Air Temperature (tas, tasmin, tasmax) Baseline + Future (Additive)

# # 1.1 Air Temperature (tas)
tas_gfdl_hist_down <- apply_bias_add_baseline("tas", model)
writeCDF(tas_gfdl_hist_down, file.path(dir_out, "tas_GFDL-ESM4_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "tas", zname = "time")
tas_gfdl_ssp245_down <- apply_bias_add_future("tas", model, "ssp245")
writeCDF(tas_gfdl_ssp245_down, file.path(dir_out, "tas_GFDL-ESM4_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "tas", zname = "time")
tas_gfdl_ssp585_down <- apply_bias_add_future("tas", model, "ssp585")
writeCDF(tas_gfdl_ssp585_down, file.path(dir_out, "tas_GFDL-ESM4_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "tas", zname = "time")

# # 1.2 Minimum Air Temperature (tasmin)
 tasmin_gfdl_hist_down <- apply_bias_add_baseline("tasmin", model)
 writeCDF(tasmin_gfdl_hist_down, file.path(dir_out, "tasmin_GFDL-ESM4_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "tasmin", zname = "time")
 tasmin_gfdl_ssp245_down <- apply_bias_add_future("tasmin", model, "ssp245")
 writeCDF(tasmin_gfdl_ssp245_down, file.path(dir_out, "tasmin_GFDL-ESM4_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "tasmin", zname = "time")
 tasmin_gfdl_ssp585_down <- apply_bias_add_future("tasmin", model, "ssp585")
 writeCDF(tasmin_gfdl_ssp585_down, file.path(dir_out, "tasmin_GFDL-ESM4_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "tasmin", zname = "time")

 # 1.3 Maximum Air Temperature (tasmax)
 tasmax_gfdl_hist_down <- apply_bias_add_baseline("tasmax", model)
 writeCDF(tasmax_gfdl_hist_down, file.path(dir_out, "tasmax_GFDL-ESM4_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "tasmax", zname = "time")
 tasmax_gfdl_ssp245_down <- apply_bias_add_future("tasmax", model, "ssp245")
 writeCDF(tasmax_gfdl_ssp245_down, file.path(dir_out, "tasmax_GFDL-ESM4_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "tasmax", zname = "time")
 tasmax_gfdl_ssp585_down <- apply_bias_add_future("tasmax", model, "ssp585")
 writeCDF(tasmax_gfdl_ssp585_down, file.path(dir_out, "tasmax_GFDL-ESM4_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "tasmax", zname = "time")

 # 2) Precipitation (pr) Baseline + Future (Logarithmic)
 pr_gfdl_hist_down <- apply_bias_log_baseline("pr", model)
 writeCDF(pr_gfdl_hist_down, file.path(dir_out, "pr_GFDL-ESM4_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "pr", zname = "time")
 pr_gfdl_ssp245_down <- apply_bias_log_future("pr", model, "ssp245")
 writeCDF(pr_gfdl_ssp245_down, file.path(dir_out, "pr_GFDL-ESM4_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "pr", zname = "time")
pr_gfdl_ssp585_down <- apply_bias_log_future("pr", model, "ssp585") # Only this one was uncommented
writeCDF(pr_gfdl_ssp585_down, file.path(dir_out, "pr_GFDL-ESM4_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "pr", zname = "time")

# # 3) Relative Humidity (hurs) Baseline + Future (Logarithmic)
hurs_gfdl_hist_down <- apply_bias_log_baseline("hurs", model) # This one was uncommented
writeCDF(hurs_gfdl_hist_down, file.path(dir_out, "hurs_GFDL-ESM4_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "hurs", zname = "time")
hurs_gfdl_ssp245_down <- apply_bias_log_future("hurs", model, "ssp245") # This one was uncommented
writeCDF(hurs_gfdl_ssp245_down, file.path(dir_out, "hurs_GFDL-ESM4_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "hurs", zname = "time")
hurs_gfdl_ssp585_down <- apply_bias_log_future("hurs", model, "ssp585") # This one was uncommented
writeCDF(hurs_gfdl_ssp585_down, file.path(dir_out, "hurs_GFDL-ESM4_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "hurs", zname = "time")

# ####### IPSL-CM6A-LR (Uncomment and set model <- "IPSL-CM6A-LR" to Run) ###############

 # 1) Air Temperature (tas, tasmin, tasmax) Baseline + Future (Additive)
 tas_ipsl_hist_down <- apply_bias_add_baseline("tas", model)
 writeCDF(tas_ipsl_hist_down, file.path(dir_out, "tas_IPSL-CM6A-LR_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "tas", zname = "time")
 tas_ipsl_ssp245_down <- apply_bias_add_future("tas", model, "ssp245")
 writeCDF(tas_ipsl_ssp245_down, file.path(dir_out, "tas_IPSL-CM6A-LR_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "tas", zname = "time")
 tas_ipsl_ssp585_down <- apply_bias_add_future("tas", model, "ssp585")
 writeCDF(tas_ipsl_ssp585_down, file.path(dir_out, "tas_IPSL-CM6A-LR_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "tas", zname = "time")

# # (TASMIN and TASMAX follow the same pattern)
# 2) Precipitation (pr) Baseline + Future (Logarithmic)
 pr_ipsl_hist_down <- apply_bias_log_baseline("pr", model)
 writeCDF(pr_ipsl_hist_down, file.path(dir_out, "pr_IPSL-CM6A-LR_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "pr", zname = "time")
 pr_ipsl_ssp245_down <- apply_bias_log_future("pr", model, "ssp245")
 writeCDF(pr_ipsl_ssp245_down, file.path(dir_out, "pr_IPSL-CM6A-LR_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "pr", zname = "time")
 pr_ipsl_ssp585_down <- apply_bias_log_future("pr", model, "ssp585")
 writeCDF(pr_ipsl_ssp585_down, file.path(dir_out, "pr_IPSL-CM6A-LR_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "pr", zname = "time")

# # 3) Relative Humidity (hurs) Baseline + Future (Logarithmic)
 hurs_ipsl_hist_down <- apply_bias_log_baseline("hurs", model)
 writeCDF(hurs_ipsl_hist_down, file.path(dir_out, "hurs_IPSL-CM6A-LR_historical_down_1979_2014.nc"), overwrite = TRUE, varname = "hurs", zname = "time")
 hurs_ipsl_ssp245_down <- apply_bias_log_future("hurs", model, "ssp245")
 writeCDF(hurs_ipsl_ssp245_down, file.path(dir_out, "hurs_IPSL-CM6A-LR_ssp245_down_2015_2100.nc"), overwrite = TRUE, varname = "hurs", zname = "time")
 hurs_ipsl_ssp585_down <- apply_bias_log_future("hurs", model, "ssp585")
 writeCDF(hurs_ipsl_ssp585_down, file.path(dir_out, "hurs_IPSL-CM6A-LR_ssp585_down_2015_2100.nc"), overwrite = TRUE, varname = "hurs", zname = "time")