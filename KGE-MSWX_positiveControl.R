# === 0. LOAD LIBRARIES ===
library(terra)      
library(tidyverse)  
library(hydroGOF)   # Contains hydrological goodness-of-fit metrics (though KGE is manually defined here)
library(lubridate)  
library(sf)         
library(purrr)      
library(ggplot2)    

# =========================================================================
# === 1. PARAMETERS AND PATH PREPARATION ===
# =========================================================================

# --- PATHS (PLACEHOLDERS) ---
# Directory containing the MSWX (Observation) data files (e.g., Relative Humidity).
dir_mswx <- "PATH_TO_MSWX_OBSERVATION_DATA/"   # <--- ADJUST THIS PATH

# Directory containing the Downscaled Climate Model (Simulation) data stacks.
dir_downscaled <- "PATH_TO_DOWNSCALED_SIMULATION_DATA/"

# Path to the administrative boundaries file for masking/cropping the analysis area.
dir_Togo <- "PATH_TO_ADMIN_BOUNDARIES/your_admin_file.gpkg"


# --- MODEL AND VARIABLE LISTS ---
models <- c(
  "GFDL-ESM4",
  "IPSL-CM6A-LR",
  "MPI-ESM1-2-HR",
  "MRI-ESM2-0",
  "UKESM1-0-LL"
)

vars <- c(
  "pr",       # Precipitation
  "tasmin",   # Min Temperature
  "tasmax",   # Max Temperature
  "tas",      # Mean Temperature
  "hurs"      # Relative Humidity
)

# Define the variable and GCM for the current run (to be iterated over)
target_variable <- vars[5]  # e.g., "hurs"
target_gcm <- models[5]     # e.g., "UKESM1-0-LL"

# --- FILE SEARCH (MSWX - Observation) ---
# Pattern example: 1979001_crop.nc (Assumes daily files are separate)
# If MSWX data is in a single stack, adjust 'mswx_files' assignment.
mswx_files <- list.files(
  path = dir_mswx,
  pattern = "_crop\\.nc$", # All MSWX files that end with '_crop.nc'
  full.names = TRUE
)

# --- FILE SEARCH (Downscaled - Simulation) ---
# Pattern example: pr_GFDL-ESM4_historical_down_1979_2014.nc (Assumes data is in a single stack)
sim_file <- list.files(
  path = dir_downscaled,
  pattern = sprintf("^%s_%s_historical_down_.*\\.nc$", target_variable, target_gcm),
  full.names = TRUE
)

if (length(mswx_files) == 0 || length(sim_file) == 0) {
  stop("No files found. Please check paths and search patterns.")
}


# =========================================================================
# === 2. KGE FUNCTION DEFINITION ===
# =========================================================================

# Vectorized function for KGE (Simulation vs. Observation)
# Applicable to time arrays (layers). 'sim' and 'obs' are time series vectors 
# for a single pixel.
kge_vec <- function(sim, obs) {
  
  # Remove NA's from valid pairs (Obs and Sim)
  valid_indices <- !is.na(sim) & !is.na(obs)
  sim_valid <- sim[valid_indices]
  obs_valid <- obs[valid_indices]
  
  N <- length(obs_valid)
  if (N < 2) return(NA_real_) # At least 2 points are needed for KGE calculation
  
  # KGE = 1 - sqrt[(r-1)^2 + (alpha-1)^2 + (beta-1)^2]
  
  # 1. Correlation (r)
  r <- cor(sim_valid, obs_valid)
  
  # 2. Variability Ratio (alpha)
  alpha <- sd(sim_valid) / sd(obs_valid)
  
  # 3. Bias Ratio (beta)
  beta <- mean(sim_valid) / mean(obs_valid)
  
  # Handle division by zero case (if obs mean is 0)
  if (is.infinite(beta) || is.nan(beta)) beta <- NA_real_ 
  
  kge_value <- 1 - sqrt((r - 1)^2 + (alpha - 1)^2 + (beta - 1)^2)
  return(kge_value)
}

# =========================================================================
# === 3. READ, ALIGN, AND CALCULATE KGE (PIXEL-BY-PIXEL) ===
# =========================================================================

# A. Read and prepare the Downscaled (Simulation) stack
r_sim <- rast(sim_file)
names(r_sim) <- time(r_sim) # Name layers by date

# B. Read and prepare the MSWX (Observation) stack
# If MSWX is in many files, create a stack (assumed to be sorted)
r_obs <- rast(mswx_files) 
names(r_obs) <- time(r_obs) # Name layers by date

# Define the aggregation function for monthly analysis
if (target_variable == "pr") {
  aggr_fun <- sum # Sum for precipitation
} else if (target_variable %in% c("tas", "tasmax", "tasmin", "hurs")) {
  aggr_fun <- mean # Mean for temperature and humidity
} else {
  stop("Variable not recognized for aggregation.")
}

# C. TEMPORAL ALIGNMENT (Monthly Aggregation)
# 1. Monthly Aggregation for Simulation (Downscaled)
sim_dates <- time(r_sim)
month_index <- format(sim_dates, "%Y-%m") 

r_sim_mensual <- terra::tapp(
  x = r_sim, 
  index = month_index, 
  fun = aggr_fun,
  filename = "", # Keep in memory
  na.rm = TRUE
)

# 2. Monthly Aggregation for Observation (MSWX)
obs_dates <- time(r_obs)
month_index_obs <- format(obs_dates, "%Y-%m")

r_obs_mensual <- terra::tapp(
  x = r_obs, 
  index = month_index_obs, 
  fun = aggr_fun,
  filename = "", 
  na.rm = TRUE
)

# C. TEMPORAL ALIGNMENT (Monthly - Clipping and Ordering)
# Clip the observation stack to match the final period of the downscaling
sim_months <- names(r_sim_mensual)
r_obs_aligned_mensual <- r_obs_mensual[[which(names(r_obs_mensual) %in% sim_months)]]

# Ensure both monthly stacks have the same layer order
r_obs_aligned_mensual <- r_obs_aligned_mensual[[order(names(r_obs_aligned_mensual))]]
r_sim_mensual <- r_sim_mensual[[order(names(r_sim_mensual))]]


# 3. SPATIAL ALIGNMENT AND MASKING
# The reference geometry is the first layer of the downscaled stack
ref_geom <- r_sim[[1]] 

# Read Admin Boundaries for final masking
togo.admin <- vect(dir_Togo)

# Resample the monthly observation stack to the original downscaling geometry
r_obs_remuestreado <- terra::resample(
  x = r_obs_aligned_mensual, 
  y = ref_geom, 
  method = "bilinear" # Bilinear interpolation for continuous variables
)

# Mutual Masking (Use the first layer of the monthly stack as a mask reference)
mascara_sim <- !is.na(r_sim_mensual[[1]])

# Apply mask based on the simulation's non-NA pixels
r_obs_final <- terra::mask(r_obs_remuestreado, mascara_sim, maskvalues = FALSE) 
r_sim_final <- terra::mask(r_sim_mensual, mascara_sim, maskvalues = FALSE)

# Crop and Mask to the Area of Interest (Togo admin boundaries)
r_obs_final <- crop(r_obs_final, ext(togo.admin))
r_obs_final <- mask(r_obs_final, togo.admin)

plot(r_obs_final[[1]]) # Check the masked observation

r_sim_final <- crop(r_sim_final, ext(togo.admin))
r_sim_final <- mask(r_sim_final, togo.admin)

plot(r_sim_final[[1]]) # Check the masked simulation

# --- KGE CALCULATION USING DATA FRAMES ---

# Get coordinates and cell IDs (Crucial for group_by operation)
r_ref <- r_sim_final[[1]]
cell_ids <- 1:ncell(r_ref) # Full list of cell IDs


# 1. Extract values from the stacks into a matrix (Nº Pixels x Nº Months)
mat_sim <- values(r_sim_final) # Simulation
mat_obs <- values(r_obs_final) # Observation

# 2. Convert matrices to "wide" data frames
df_sim <- as.data.frame(mat_sim)
df_obs <- as.data.frame(mat_obs)


# 3. Add Pixel IDs and transform to "long" format
df_largo_sim <- df_sim %>%
  mutate(CellID = cell_ids, .before = 1) %>%
  pivot_longer(
    cols = -CellID, 
    names_to = "MonthYear", 
    values_to = "Sim"
  )

df_largo_obs <- df_obs %>%
  mutate(CellID = cell_ids, .before = 1) %>%
  pivot_longer(
    cols = -CellID, 
    names_to = "MonthYear", 
    values_to = "Obs"
  )

# 4. Join Data Frames by Pixel and Month
df_kge_calc <- inner_join(df_largo_sim, df_largo_obs, by = c("CellID", "MonthYear")) %>%
  # Remove NA's (only valid pairs are needed for KGE calculation)
  na.omit()


# 5. KGE Calculation grouped by CellID
df_kge_resultado <- df_kge_calc %>%
  group_by(CellID) %>%
  summarise(
    # Apply the KGE function to the time series of each pixel
    KGE_Value = kge_vec(sim = Sim, obs = Obs),
    .groups = "drop"
  )

print(df_kge_resultado)


# 6. Create the Final KGE SpatRaster

# Create a vector of KGE values, ordered by CellID
kge_vector <- rep(NA_real_, ncell(r_ref)) # Initialize with NA's for all pixels
kge_vector[df_kge_resultado$CellID] <- df_kge_resultado$KGE_Value # Fill calculated pixels

# Create the new raster from a reference layer
kge_map_df <- r_ref
values(kge_map_df) <- kge_vector

# Assign Final Name
names(kge_map_df) <- paste0("KGE_", target_variable, "_", target_gcm, "_DF")
print(kge_map_df)

# Visualization (for the current run)
plot(kge_map_df, main = "KGE Pixel-by-Pixel (Data Frame Method)")


# --- POST-CALCULATION STEPS (for iteration/saving results) ---

comb.name <- paste0("KGE_", target_variable, "_", target_gcm)

# NOTE: The list 'list.KGE.results' must be initialized before running the loop.
# list.KGE.results <- list()
# This block is intended to be run *within* a loop over all variables/GCMs.
# list.KGE.results <- c(list.KGE.results, list(df_kge_resultado))
# names(list.KGE.results)[[length(list.KGE.results)]] <- comb.name

# After evaluating every model-var, bind list rows
# df.kge.results.final <- bind_rows(list.KGE.results, .id = "var_GCM")


# =========================================================================
# === Prepare Data for Plotting (Assumes all models/vars have been run) ===
# =========================================================================

# The following assumes df.kge.results.final has been created by binding all results.

# 1. Apply strsplit to the combined column
# split_list <- strsplit(df.kge.results.final$var_GCM, split = "_")

# 2. Extract elements to create separate columns for Variable and GCM
# df.kge.results.final <- df.kge.results.final %>%
#   mutate(
#     Variable = map_chr(split_list, ~ .x[2]),
#     GCM = map_chr(split_list, ~ .x[3])
#   )

# print(head(df.kge.results.final))

# =========================================================================
# === 4. CONVERSION AND SAVING FUNCTION ===
# =========================================================================

# Function to convert the KGE_Value back to a raster and prepare for saving/plotting
crear_y_graficar_kge_map <- function(df_grupo, r_ref_geom) {
  
  var_name <- unique(df_grupo$Variable)
  gcm_name <- unique(df_grupo$GCM)
  
  kge_vector <- rep(NA_real_, ncell(r_ref_geom)) 
  kge_vector[df_grupo$CellID] <- df_grupo$KGE_Value
  
  kge_map_final <- r_ref_geom
  values(kge_map_final) <- kge_vector
  
  # The function returns the SpatRaster
  return(kge_map_final) 
}

# 1. Group the final results data frame by Variable and GCM
# df_grupos <- df.kge.results.final %>%
#   mutate(CellID = as.numeric(CellID)) %>%
#   group_by(Variable, GCM)

# 2. Generate the list of SpatRasters
# lista_de_mapas <- df_grupos %>%
#   group_map(~ crear_y_graficar_kge_map(., r_ref_geom = r_ref)) %>%
#   
#   # 3. Assign names to the list
#   set_names(
#     nm = df_grupos %>% 
#       group_keys() %>% 
#       unite(col = "MapName", Variable, GCM, sep = "_") %>% 
#       pull(MapName)
#   )

# print(names(lista_de_mapas))

# kge.rasters <- rast(lista_de_mapas)

# --- SAVE RESULTS (PLACEHOLDERS) ---
# setwd("PATH_TO_KGE_RESULTS_DIRECTORY/")
# saveRDS(kge.rasters, "kge_rasters_resultados.RDS")
# saveRDS(df.kge.results.final, "kge_table_results.RDS")