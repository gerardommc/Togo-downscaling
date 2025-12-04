# =========================================================================
# DAILY BIAS (DOY) CALCULATION SCRIPT
# This script computes the Day Of Year (DOY) climatological bias 
# between a coarse GCM (Attrici/NEX-GDDP) and a reference dataset (GloH2O).
# The calculated bias is used for subsequent bias correction and downscaling.
# =========================================================================

library(terra)
library(sf)
library(fields) # For Thin Plate Spline (TPS) in later steps (or bias calculation itself)
library(lubridate)
library(ncdf4) # Useful for low-level NetCDF operations, though terra handles most.

##### CONFIGURATION (ADJUST THESE) ---------------------------------------

# Variables to process
vars <- c("pr", "tasmin", "tasmax", "tas", "hurs")

# Select the variable to process (Currently set to the 3rd element: "tasmax")
var <- vars[3] # "pr" | "tasmin" | "tasmax" | "tas" | "hurs"

# Baseline reference period for bias calculation
period_ref <- c(1979, 2014)

# Geometries and Grid Paths
aoi_path <- "PATH_TO_ADMIN_BOUNDARIES/gadm_TGO_buffer.gpkg"
grid_path <- "PATH_TO_ADMIN_BOUNDARIES/grid10km.nc" # Your target 10km grid template (with values in AOI)

# Reference Observation Data (GloH2O)
# Directory containing GloH2O NetCDF files (1 .nc per day)
gloh_dir <- file.path("PATH_TO_REFERENCE_DATA/GloH2O/tmax_recortados/") 

# Model Data (Attrici/NEX-GDDP)
# Directory listing model files (coarse resolution, multi-year)
attrici.dirs <- list.dirs("PATH_TO_MODEL_DATA/attrici/")
attrici_files <- list.files(attrici.dirs, pattern="\\.nc$", full.names = TRUE)

# Filter model files for the selected variable (var)
hist.var.attrici.files <- attrici_files[grepl(var, attrici_files)]

# Output file definition
out.dir <- "PATH_TO_BIAS_OUTPUT_DIRECTORY/attrici/"
# Define the output NetCDF filename using the selected variable
bias_out_nc <- paste(out.dir, "attrici_", var, "_", period_ref[1], "_", period_ref[2], "_Bias_DOY.nc", sep = "")

# Rolling window size for DOY climatology (e.g., ± 15 days)
roll_k <- 15 

# =========================
# UTILITY FUNCTIONS
# =========================

# Function to safely extract time information from a NetCDF path using terra
extract_terra_date <- function(nc_path) {
  r <- tryCatch({ rast(nc_path) }, error = function(e) {
    warning(paste("Could not read file with terra:", nc_path)); return(NULL) })
  if (is.null(r)) { return(as.Date(NA)) }
  return(as.Date(time(r)))
}

# Function to find indices corresponding to a rolling Day Of Year window (±k)
doy_window_idx <- function(dates_vec, k = 15) {
  d <- yday(dates_vec) # 1..366
  lapply(1:366, function(D) {
    # Calculate circular difference (handles year boundaries)
    dif <- ((d - D + 366 + 183) %% 366) - 183
    # Return indices where the absolute difference is within k days
    which(abs(dif) <= k)
  })
}

# =========================
# AOI / GRID SETUP
# =========================

# Load and transform AOI geometry
aoi_sf <- st_read(aoi_path, quiet=TRUE) |> st_transform(4326)
aoi <- vect(aoi_sf)

# Load the 10km grid template
grid01_mask <- rast(grid_path); crs(grid01_mask) <- "EPSG:4326"
if (!hasValues(grid01_mask)) grid01_mask[] <- 1L

# =========================
# 1. LOAD MODEL DATA (Attrici/NEX multi-year) → 1979–2014
# =========================
mod_day <- rast(hist.var.attrici.files); crs(mod_day) <- "EPSG:4326"

dates_m <- time(mod_day)

# Check if time metadata is present (common issue with TIFs/some NetCDFs)
if (is.null(dates_m) || all(is.na(dates_m))) {
  # --- Fallback for data without internal time metadata ---
  
  # **REQUIRED ADJUSTMENT: Set 'origin_date' to the TRUE start date of your time series**
  # Example: If your first file starts Jan 1, 1950.
  origin_date <- as.Date("1950-01-01") 
  
  num_layers <- nlyr(mod_day)
  
  # Calculate dates: origin + day sequence (0-indexed)
  dates_m <- origin_date + seq_len(num_layers) - 1
  
} else {
  # Ensure dates are explicitly of Date class
  dates_m <- as.Date(dates_m)
}

# Select layers corresponding to the reference period (1979-2014)
selm <- which(year(dates_m) >= period_ref[1] & year(dates_m) <= period_ref[2])
mod_day <- mod_day[[selm]]; dates_m <- dates_m[selm]

# Crop model data to the extent of the 10km grid template
mod_day <- crop(mod_day, ext(grid01_mask) )

# NASA UNITS CONVERSION (if data is in Kelvin or strange units)
if (var == "pr") {
  # pr: kg m^-2 s^-1 → mm/day
  mod_day <- mod_day * 86400
} else if (var %in% c("tas", "tasmax", "tasmin")) {
  # temps: Kelvin → °C
  mod_day <- mod_day - 273.15
}
# =========================
# 2. LOAD OBSERVATION DATA (GloH2O daily files) → 1979–2014
# =========================
gloh_paths <- list.files(gloh_dir, pattern="\\.nc$", full.names=TRUE)

# Extract dates from each GloH2O NetCDF path/metadata
dates_g <- map_chr(gloh_paths, extract_terra_date) |> as.Date()

# Select GloH2O files corresponding to the reference period
selg <- which(year(dates_g) >= period_ref[1] & year(dates_g) <= period_ref[2])
gloh_paths <- gloh_paths[selg]; dates_g <- dates_g[selg]

# Ensure adequate overlap/data existence
stopifnot(length(dates_g) > 0, length(dates_m) > 0)

# =========================
# 3. DOY CLIMATOLOGIES (±k) ON MODEL GRID
# =========================

# 3.1 DOY Climatology for Model Data (Attrici/NEX)
idx_doy_m <- doy_window_idx(dates_m, k = roll_k)
mod_tpl <- mod_day[[1]] # Template layer for model resolution
mod_doy <- rast(mod_tpl, nlyr=366) # Output DOY climatology raster (366 layers)

# Compute mean for each DOY rolling window
for (D in 1:366) {
  ii <- idx_doy_m[[D]]
  mod_doy[[D]] <- if (length(ii) > 0) mean(mod_day[[ii]], na.rm=TRUE) else NA * mod_tpl
}

# 3.2 DOY Climatology for Observation Data (GloH2O) on Model Grid
gloh_doy_on_mod <- rast(mod_tpl, nlyr=366); gloh_doy_on_mod[] <- NA
sum_coarse <- lapply(1:366, function(i) mod_tpl * 0) # List to hold cumulative sums
cnt_coarse <- integer(366) # Vector to hold counts

# Aggregate GloH2O data by DOY
for (i in seq_along(gloh_paths)) {
  r <- rast(gloh_paths[i]); crs(r) <- "EPSG:4326"
  r <- crop(r, ext(mod_tpl))
  
  # GloH2O Units: Assume units are already standard (mm/day or °C). 
  # Apply conversion here if needed.
  
  D <- yday(dates_g[i])
  # Resample GloH2O (fine) to the coarse Model grid resolution
  r_on_mod <- resample(r, mod_tpl, method="bilinear")
  
  # Accumulate sum and count for the respective DOY
  sum_coarse[[D]] <- cover(sum_coarse[[D]] + r_on_mod, sum_coarse[[D]])
  cnt_coarse[D]   <- cnt_coarse[D] + 1
}

# Calculate the mean for each DOY
for (D in 1:366) {
  gloh_doy_on_mod[[D]] <- if (cnt_coarse[D] > 0) sum_coarse[[D]] / cnt_coarse[D] else NA * mod_tpl
}
rm(sum_coarse, cnt_coarse)

# =========================
# 4. CALCULATE DOY BIAS (Additive or Logarithmic)
# =========================

# Epsilon value for log transform (avoids log(0))
if (var == "pr") {
  eps <- 0.1     # in mm/day
} else if (var == "hurs") {
  eps <- 0.01    # for humidity (adjust if using 0–1 range)
} else {
  eps <- 0       # for tas* (no log required)
}

# Compute Bias
if (var %in% c("pr", "hurs")) {
  # Logarithmic Bias (for variables where negative values are physically impossible/problematic)
  # Bias = log(Model/Obs) = log(Model) - log(Obs)
  mod_doy_log <- log(mod_doy + eps)
  gloh_doy_log <- log(gloh_doy_on_mod + eps)
  bias_doy_on_mod <- mod_doy_log - gloh_doy_log # delta_log
} else {
  # Additive Bias (for temperature variables)
  # Bias = Model - Obs
  bias_doy_on_mod <- mod_doy - gloh_doy_on_mod
}

# =========================
# 5. SAVE DOY BIAS
# =========================

# The result is a SpatRaster stack with 366 layers (DOY 1 to 366)
writeCDF(bias_doy_on_mod, bias_out_nc, overwrite=TRUE)
message("DOY Bias calculation complete and saved to: ", bias_out_nc)