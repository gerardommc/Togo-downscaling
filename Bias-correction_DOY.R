### Bias NASA dataset

#1) 01_bias_doy.R — Daily bias calculation (DOY)

# 01_bias_doy.R
# Daily bias (DOY ±15) NASANEX vs GloH2O in the model grid
# Output: NetCDF file with 366 DOY bias layers (factor for pr; delta for tas)
# This script was developed by Sandra Milena Castaño Quintero, https://github.com/Salamanadra


pkg <- c("terra", "sf", "fields", "lubridate", "ncdf4")
sapply(pkg, \(p) if (!requireNamespace(p, quietly = TRUE)) install.packages(p))
sapply(pkg, library, character.only = TRUE)

# ========================================================================================================
# PARAMETERS
# ========================================================================================================
# Climmatic variables: "pr" - precipitation | "tasmin" - air. temp. min. | "tasmax"  - air temp. max |"tas" - air temperature | "hurs" - relative humidity
# 5 Global Climate Models - GCMs
models <- c(
  "GFDL-ESM4",
  "IPSL-CM6A-LR",
  "MPI-ESM1-2-HR",
  "MRI-ESM2-0",
  "UKESM1-0-LL")

vars<-c(
  "pr",
  "tasmin",
  "tasmax",
  "tas",
  "hurs")

# Directory
period_ref  <- c(1979, 2014)
aoi_path    <- "E:/RS/gadm_TGO_buffer.gpkg"
grid_path   <- "E:/RS/grid10km.nc"         #grid 0.1° (with values in AOI)
# Directory Historical baseline 
gloh_dir    <- file.path("E:/RS/GloH2O/Historic_croppeds/pr_cropped/") # GloH2O: 1 .nc x day
# Directory - GCMS / Future scenario
nasa.dirs   <- list.dirs("E:/RS/NASA-NEX/NASA/")
nasa_files  <- list.files(nasa.dirs, pattern="\\.tif$",full.names = T)
# Directory - Counterfactual scenario
attrici.dirs<-list.dirs("D:/RS/COUNTERCLIM/cropped/")  #For Counterfactual Model, we use Atricci  dataset
attrici_files <- list.files(attrici.dirs, pattern="\\.nc$",full.names = T)

# vars value [1,2,3,4,5] will be changed to the one we are working on from the line 25-30;  e.g.: "tasmax" = [3]
# models value [1,2,3,4,5] will be changed to the one we are working on from the line 18-23; e.g.: "MRI-ESM2-0" = [4]
hist.var.gcm.files<-nasa_files[grepl(vars[1],nasa_files)&grepl(models[5],nasa_files)&grepl("historical",nasa_files)]
# For counterfactual model, we used atricci dataset, aviable in "atricci_files" object, only change var value (line 25-30).
hist.var.attrici.files<-attrici_files[grepl(vars[2],attrici_files)]

files_to_process_gcms <- hist.var.gcm.files 
files_to_process_counterfactual <- hist.var.attrici.files

out.dir<-"MODELS/bias_DOY/" 

#change file title title to the variable and model we are working with: 
#      e.g.: "/bias_%dataset%_%baseline%_%var%_%model%_%period%_Bias_DOY.nc" 
bias_out_nc<-paste(out.dir,"/bias_NASA_baseline_pr_UKESM1-0-LL_1979_2014_Bias_DOY.nc",sep ="")

roll_k <- 15   # window DOY ±k days

# ========================================================================================================
# UTILITIES
# ========================================================================================================
extract_terra_date <- function(nc_path) {
  r <- tryCatch({ rast(nc_path) }, error = function(e) {
    warning(paste("Could not read file with terra:", nc_path)); return(NULL) })
  if (is.null(r)) { return(as.Date(NA)) }
  return(as.Date(time(r)))
}

doy_window_idx <- function(dates_vec, k = 15) {
  d <- yday(dates_vec)              # 1..366
  lapply(1:366, function(D) {
    dif <- ((d - D + 366 + 183) %% 366) - 183
    which(abs(dif) <= k)
  })
}

# ============================================================================================================
# AOI / GRID
# ============================================================================================================
aoi_sf <- st_read(aoi_path, quiet=TRUE) |> st_transform(4326)
aoi    <- vect(aoi_sf)

grid01_mask <- rast(grid_path); crs(grid01_mask) <- "EPSG:4326"
if (!hasValues(grid01_mask)) grid01_mask[] <- 1L

# ============================================================================================================
# LOAD NASA (multi-year) → 1979–2014
# ============================================================================================================
mod_day <- rast(hist.var.gcm.files); crs(mod_day) <- "EPSG:4326"
# If 0..360: mod_day <- wrap(mod_day, "lon", xmin=-180)

dates_m <- time(mod_day)

# Verify if "mod_day" (raster/terra object)  has time information (comon for TIF is NULL)
if (is.null(dates_m) || all(is.na(dates_m))) {
  # --- Fallback for TIF without internal time metadata 
  
  # 1. Assume the NEX-GDDP convention: Origin and consecutive days.
  #   Requires you know the start date of the first file or series.
  
  # Example: If you are reading a series that begins on January 1, 1950:
  # **ADJUST 'origin_date' to the REAL start date of your TIF series.
  origin_date <- as.Date("1950-01-01") 
  
  # Determine the number of layers (days) in the object 'mod_day' (line 84).
  num_layers <- nlyr(mod_day)
  
  # Calculates dates: origin + sequence of days (0-indexed).
  dates_m <- origin_date + seq_len(num_layers) - 1
  
} else {
  # If 'mod_day' object had dates (e.g., was created from NetCDF), make sure those are of the type "Date".
  dates_m <- as.Date(dates_m)
}
selm    <- which(year(dates_m) >= period_ref[1] & year(dates_m) <= period_ref[2])
mod_day <- mod_day[[selm]]; dates_m <- dates_m[selm]
mod_day <- crop(mod_day, ext(grid01_mask) )

####### UNITS in NASA dataset ###### 
# When the processed variable is precipitation, the values are transformed from kg m⁻² s⁻¹ to mm day⁻¹ by multiplying by 86,400;
# if it is a temperature variable (tas, tasmax, or tasmin), the units are converted from kelvin to degrees Celsius by subtracting 273.15.

if (var == "pr") {
  mod_day <- mod_day * 86400
} else if (var %in% c("tas", "tasmax", "tasmin")) {
  mod_day <- mod_day - 273.15
}

#  ============================================================================================================
# LOAD GloH2O dataset dayli (one archive x day) → 1979–2014
#  ============================================================================================================
gloh_paths <- list.files(gloh_dir, pattern="\\.nc$", full.names=TRUE)
dates_g    <- extract_terra_date(gloh_paths)
selg       <- which(year(dates_g) >= period_ref[1] & year(dates_g) <= period_ref[2])
gloh_paths <- gloh_paths[selg]; dates_g <- dates_g[selg]

# Ensure adequate overlap
stopifnot(length(dates_g) > 0, length(dates_m) > 0)

# ===============================================================================================================
# DOY CLIMATOLOGY (±k) in grid's Model.
# ===============================================================================================================
# 1) DOY NASA (coarse)
idx_doy_m <- doy_window_idx(dates_m, k = roll_k)
mod_tpl <- mod_day[[1]]
mod_doy <- rast(mod_tpl, nlyr=366)
for (D in 1:366) {
  ii <- idx_doy_m[[D]]
  mod_doy[[D]] <- if (length(ii)>0) mean(mod_day[[ii]], na.rm=TRUE) else NA*mod_tpl
}

# 2) DOY GloH2O (coarse) – average by DOY reading day by day.
gloh_doy_on_mod <- rast(mod_tpl, nlyr=366); gloh_doy_on_mod[] <- NA
sum_coarse <- lapply(1:366, function(i) mod_tpl*0)
cnt_coarse <- integer(366)

for (i in seq_along(gloh_paths)) {
  r <- rast(gloh_paths[i]); crs(r) <- "EPSG:4326"
  r <- crop(r, ext(mod_tpl))
  #GloH2O units: pr in mm/day; temps in °C. If not, convert here.
  D <- yday(dates_g[i])
  r_on_mod <- resample(r, mod_tpl, method="bilinear")
  sum_coarse[[D]] <- cover(sum_coarse[[D]] + r_on_mod, sum_coarse[[D]])
  cnt_coarse[D]   <- cnt_coarse[D] + 1
}
for (D in 1:366) {
  gloh_doy_on_mod[[D]] <- if (cnt_coarse[D]>0) sum_coarse[[D]]/cnt_coarse[D] else NA*mod_tpl
}
rm(sum_coarse, cnt_coarse)

# Constant to avoid log(0). Adjust according "pr" o "hurs" variable.
if (var == "pr") {
  eps <- 0.1      # mm/day
} else if (var == "hurs") {
  eps <- 0.01     # if var "hurs" it's in 0–100; then use 0.0001 if it's en 0–1
} else {
  eps <- 0        # for var "tas*" it is not necessary,
}

# =====================================================================================================
# BIAS DOY IN MODEL GRID (logarithmic difference for "pr" and "hurs" variables).
# =====================================================================================================
if (var %in% c("pr", "hurs")) {
  # Bias as subtraction of logarithms
  # (for hurs, prevents corrections from giving negative values when applying bias).
  mod_doy_log  <- log(mod_doy + eps)
  gloh_doy_log <- log(gloh_doy_on_mod + eps)
  bias_doy_on_mod <- mod_doy_log - gloh_doy_log   # delta_log
} else {
  # Tas*: additive bias (delta).
  bias_doy_on_mod <- mod_doy - gloh_doy_on_mod
}
# ========================================================================================================
# SAVE BIAS DOY
# ========================================================================================================
bias_out_nc
writeCDF(bias_doy_on_mod, bias_out_nc, overwrite=TRUE)
message("Ready: ", bias_out_nc)