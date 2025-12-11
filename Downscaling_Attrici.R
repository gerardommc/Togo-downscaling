# =========================================================================
# ATTRICI (Counterclim) Downscaling and Bias Correction
# Applies DOY bias (GloH2O vs ATTRICI) and interpolates to 10 km grid.
# Target period: Historical Baseline (1979–2014)
# This script was developed by Sandra Milena Castaño Quintero, https://github.com/Salamanadra
# =========================================================================

library(terra)
library(lubridate)
library(dplyr)
library(fields) # for Tps (Thin Plate Spline)
library(stringr) # for string manipulation

##### DIRECTORY PATHS (ADJUST THESE) ---------------------------------------

# ATTRICI (counterclim) historical TIFs/NetCDFs (coarse resolution)
# These files should cover the 1979-2014 period.
dir_nex <- "PATH_TO_ATTRICI_CROPPED_INPUT/"

# NetCDF files with pre-calculated DOY bias (GloH2O vs ATTRICI coarse)
# E.g.: attrici_tas_1979_2014_Bias_DOY.nc
dir_bias <- "PATH_TO_DOY_BIAS_ATTRICI_INPUT/"

# Custom 10 km grid (target fine grid template)
grid10km_path <- "PATH_TO_CUSTOM_GRID_TEMPLATE/grid10km.nc"

# Output directory for corrected and downscaled NetCDF files
dir_out <- "PATH_TO_DOWNSCALED_ATTRICI_OUTPUT/"
if (!dir.exists(dir_out)) dir.create(dir_out)

# Model identifier (fixed for ATTRICI counterfactual run)
model <- "counterclim"

# Baseline period for processing (should match the bias calculation period)
baseline_start <- as.Date("1979-01-01")
baseline_end <- as.Date("2014-12-31")

# Epsilon value for logarithmic correction (adjust if necessary)
eps_pr <- 0.1  # mm/day (Precipitation)
eps_hurs <- 0.01 # Humidity (adjust if input is 0–1 instead of 0–100)

# Variables to process (tasmax is currently selected in the execution loop)
variables_to_process <- c("tas", "tasmin", "tasmax", "pr", "hurs")

##### LOAD 10 KM GRID TEMPLATE -------------------------------------

grid10km <- rast(grid10km_path)
crs(grid10km) <- "EPSG:4326"

if (!hasValues(grid10km)) grid10km[] <- 1L

##### AUXILIARY FUNCTIONS (TPS and DOY) -----------------------------

doy_vec <- function(dates) yday(dates)

# Thin Plate Spline (TPS) interpolation from coarse layer to fine grid
tps_to_fine <- function(bias_layer, template_fine) {
  vals <- values(bias_layer, mat = FALSE)
  ok <- !is.na(vals)
  if (!any(ok)) {
    warning("All cells in bias_layer are NA; returning NA in the fine layer.")
    r_out <- template_fine
    values(r_out) <- NA
    return(r_out)
  }
  cells_ok <- which(ok)
  pts_ok <- terra::xyFromCell(bias_layer, cells_ok)
  vals_ok <- vals[ok]
  
  # Fit TPS model
  tps_fit <- fields::Tps(pts_ok, vals_ok)
  
  # Predict onto the fine grid coordinates
  fine_cells <- seq_len(ncell(template_fine))
  fine_xy <- terra::xyFromCell(template_fine, fine_cells)
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
# Load ATTRICI Baseline (Historical)
# ---------------------------------
load_attrici_baseline <- function(var) {
  
  # 1. SEARCH FOR FILES
  # Pattern ensures exclusive search for the desired variable (e.g., 'tas', not 'tasmin')
  new_pattern <- sprintf("^gswp3v109-w5e5_counterclim_%s_global_daily_.*_Crop\\.nc$", var)
  
  files_nex <- list.files(
    dir_nex,
    pattern = new_pattern,
    full.names = TRUE
  )
  
  if (length(files_nex) == 0) {
    stop(paste("No ATTRICI files found for variable:", var))
  }
  
  # 2. CHRONOLOGICAL ORDERING
  # Extract the start year (AAAA) from the AAAA_BBBB range in the filename for sorting
  year_starts <- str_sub(str_extract(basename(files_nex), "\\d{4}_\\d{4}"), 1, 4)
  files_nex <- files_nex[order(as.numeric(year_starts))]
  
  message(paste("Files found and ordered for", var, ":", length(files_nex)))
  
  # 3. LOAD ALL CHUNKS INTO A SINGLE STACK
  r_all <- rast(files_nex)
  
  
  # 4. UNITS CONVERSION
  if (var == "pr") {
    # pr: kg m^-2 s^-1 → mm/day
    r_all <- r_all * 86400
  } else if (var %in% c("tas", "tasmin", "tasmax")) {
    # temps: Kelvin → °C
    r_all <- r_all - 273.15
  }
  
  # 5. DATE ASSIGNMENT (Crucial step)
  
  # Extract start year from the first ordered file (e.g., 1911 from 1911_1920)
  start_year_loaded <- as.numeric(str_sub(str_extract(basename(files_nex[1]), "\\d{4}_\\d{4}"), 1, 4))
  
  # Define the full sequence range based on ATTRICI coverage
  dates_start <- as.Date(paste0(start_year_loaded, "-01-01"))
  dates_end <- as.Date("2019-12-31") # ATTRICI typically runs up to 2019
  
  dates_full_range <- seq(dates_start, dates_end, by = "day")
  
  # Trim the date sequence to match the actual number of loaded layers
  dates_to_assign <- dates_full_range[seq_len(nlyr(r_all))]
  
  if (anyNA(dates_to_assign)) {
    stop("Internal date sequence error. Number of layers does not match expected dates.")
  }
  
  time(r_all) <- dates_to_assign
  
  # 6. SELECT THE BASELINE PERIOD (1979-2014)
  sel <- time(r_all) >= baseline_start & time(r_all) <= baseline_end
  
  if (sum(sel) == 0) {
    stop("Error: No layers found in the 1979-2014 date range after selection.")
  }
  
  return(r_all[[sel]])
}

# ---------------------------------
# Load DOY Bias (coarse) and interpolate to 10 km grid
# ---------------------------------
load_bias_doy_fine <- function(var, template_fine) {
  # The bias file name must match the output of the previous bias calculation script
  bias_file <- file.path(
    dir_bias,
    sprintf("attrici_%s_1979_2014_Bias_DOY.nc", var)
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

##### APPLY ADDITIVE BIAS CORRECTION (tas, tasmin, tasmax) --------
# Corrected = ATTRICI_fine - Bias_DOY_fine

apply_bias_add_baseline <- function(var) {
  # Load ATTRICI, Resample to fine grid, Load Bias DOY
  r_attrici <- load_attrici_baseline(var)
  dates_b <- time(r_attrici)
  template_fine <- load_template_fine(var) # grid10km
  r_attrici_fine <- resample(r_attrici, template_fine, method = "bilinear")
  time(r_attrici_fine) <- dates_b
  doy_b <- yday(dates_b)
  nb_b <- nlyr(r_attrici_fine)
  bias_fine <- load_bias_doy_fine(var, template_fine)
  
  out_list <- vector("list", nb_b)
  for (i in 1:nb_b) {
    d <- doy_b[i]
    bias_layer <- bias_fine[[d]]
    attrici_layer <- r_attrici_fine[[i]]
    
    # ADDITIVE CORRECTION: ATTRICI_corr = ATTRICI - (ATTRICI - OBS)
    out_list[[i]] <- attrici_layer - bias_layer
  }
  
  r_corr <- rast(out_list)
  time(r_corr) <- dates_b
  
  # --- CRUCIAL STEP: APPLY FINAL MASK ---
  # Ensures the corrected data is only present within the 10km grid AOI
  r_masked <- terra::mask(r_corr, template_fine)
  
  return(r_masked)
}

##### APPLY LOGARITHMIC BIAS CORRECTION (pr, hurs) ----------------------
# Corrected = exp( log(ATTRICI + eps) - log(Bias)) - eps

apply_bias_log_baseline <- function(var) {
  # Load ATTRICI, Resample to fine grid, Load Bias DOY
  eps <- if (var == "pr") eps_pr else if (var == "hurs") eps_hurs else 0
  r_attrici <- load_attrici_baseline(var)
  dates_b <- time(r_attrici)
  template_fine <- load_template_fine(var) # grid10km
  r_attrici_fine <- resample(r_attrici, template_fine, method = "bilinear")
  time(r_attrici_fine) <- dates_b
  doy_b <- yday(dates_b)
  nb_b <- nlyr(r_attrici_fine)
  log_bias_fine <- load_bias_doy_fine(var, template_fine)
  
  out_list <- vector("list", nb_b)
  for (i in 1:nb_b) {
    d <- doy_b[i]
    attrici_layer <- r_attrici_fine[[i]]
    
    # Clip negative values before log transform
    vals <- values(attrici_layer); vals[vals < 0] <- 0; values(attrici_layer) <- vals
    
    log_attrici <- log(attrici_layer + eps)
    log_bias_day <- log_bias_fine[[d]]
    
    # LOGARITHMIC CORRECTION: log(OBS) = log(ATTRICI) - log(Bias)
    log_corr <- log_attrici - log_bias_day
    var_corr <- exp(log_corr) - eps # Inverse transform
    
    # Clip corrected negative values
    vals_corr <- values(var_corr); vals_corr[vals_corr < 0] <- 0; values(var_corr) <- vals_corr
    
    out_list[[i]] <- var_corr
  }
  
  r_corr <- rast(out_list)
  time(r_corr) <- dates_b
  
  # --- CRUCIAL STEP: APPLY FINAL MASK ---
  r_masked <- terra::mask(r_corr, template_fine)
  
  return(r_masked)
}


##### MAIN EXECUTION LOOP -----------------------------
message("Starting Bias Correction and Downscaling for ATTRICI (1979-2014)...")

# Note: The variables_to_process variable was reset to only include "tasmax" in the original script.
# We will use the full list defined at the top for maximum utility, but you can uncomment the line below:
# variables_to_process <- c("tasmax") 

for (var in variables_to_process) {
  
  message(paste("\n--- Processing variable:", var, "---"))
  
  if (var %in% c("tas", "tasmin", "tasmax")) {
    
    # 1. Apply Additive Correction (Temperatures)
    r_corrected <- tryCatch({
      apply_bias_add_baseline(var)
    }, error = function(e) {
      message(paste("❌ ERROR processing", var, ":", e$message)); return(NULL)
    })
    
  } else if (var %in% c("pr", "hurs")) {
    
    # 1. Apply Logarithmic Correction (Precipitation and Humidity)
    r_corrected <- tryCatch({
      apply_bias_log_baseline(var)
    }, error = function(e) {
      message(paste("❌ ERROR processing", var, ":", e$message)); return(NULL)
    })
    
  } else {
    message(paste("⚠️ Variable", var, "omitted (not supported)."))
    next
  }
  
  # 2. Save corrected SpatRaster to NetCDF
  if (!is.null(r_corrected)) {
    output_filename <- file.path(
      dir_out, 
      sprintf("%s_%s_down_1979_2014.nc", var, model)
    )
    
    writeCDF(
      r_corrected,
      output_filename,
      overwrite = TRUE,
      varname = var,
      zname = "time"
    )
    
    message(paste("✅ Saved:", basename(output_filename)))
  }
}

message("\n--- ATTRICI Downscaling and Bias Correction Process Finished ---")