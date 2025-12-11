# ==============================================================================
# SCRIPT NAME: Script_Togo_GCM_Validation_Monthly_KGE.R
# DESCRIPTION: 
#   1. Processes downscaled files (Historical 1979-2014).
#   2. Cleans and aggregates Togo station data (precip, temp) to monthly means.
#   3. Validates GCM performance against stations using KGE (Kling-Gupta Efficiency).
#   4. Estimates uncertainty using Bootstrap resampling (95% CI).
#
# INPUTS: 
#   - Downscaled NetCDF files (NASA NEX-GDDP-CMIP6)
#   - Daily Station CSV files (Togo meteo network)
#
# DATE:   December 2025
# This script was developed by Sandra Milena Castaño Quintero, https://github.com/Salamanadra
# ==============================================================================

##### create a database with the monthly averages for each year #####

# Load required packages
library(terra)      # For handling spatial raster data (NetCDF)
library(tidyverse)  # For data manipulation (dplyr, tibble) and clean coding
library(lubridate)  # For easy date and time operations
library(ggplot2)    
library(scales)     
library(viridis)    
library(RColorBrewer)


# --- Define Global Parameters ---
# **CONFIGURATION REQUIRED:** Set the folder where your NetCDF files are located
PATH_TO_DOWNSCALED_NC <- "PATH_TO_downscaled_nasa_2/" 

# List of climate variables to process
vars<-c("pr",     
        "tas",    
        "tasmin", 
        "tasmax", 
        "hurs")    

# List of Global Climate Models (GCMs)
models<-c("GFDL-ESM4",
          "IPSL-CM6A-LR",
          "MPI-ESM1-2-HR",
          "MRI-ESM2-0",
          "UKESM1-0-LL")

# Pattern to find all relevant NetCDF files (*_historical_down_*.nc)
all_nc_files <- list.files(
  path = PATH_TO_DOWNSCALED_NC, 
  pattern = "_historical_down_.*\\.nc$",
  full.names = TRUE, 
  recursive = FALSE # Files are directly in the specified folder
)


# --- Core Function to Process a Single NetCDF File ---
# This function calculates the spatial mean and aggregates to monthly average.
process_nc_terra_month <- function(nc_file) {
  
  # 1. Extract Metadata (Variable and GCM) from the filename
  base_name <- basename(nc_file)
  parts <- unlist(strsplit(base_name, "_"))
  
  variable <- parts[1]
  gcm_model <- parts[2]
  
  # 2. Open the NetCDF file as a SpatRaster (stack of daily rasters)
  r <- rast(nc_file)
  
  # 3. Calculate the Daily Spatial Average (mean across all grid cells for each day)
  df_daily_spatial_mean <- global(r, fun = "mean", na.rm = TRUE)
  
  # 4. Extract the time/date information
  dates <- time(r)
  
  # 5. Create the Daily Data Frame with metadata
  df_daily_result <- data.frame(
    Time = dates,
    Mean = df_daily_spatial_mean[, 1], # Spatial daily average value
    Variable = variable,
    GCM = gcm_model
  )
  
  # 6. Grouping and Calculation of Monthly Average
  df_monthly_result <- df_daily_result %>%
    mutate(
      Year = year(Time),
      Month = month(Time, label = FALSE) # Numeric month (1-12)
    ) %>%
    # Group by the 4 identifier columns (Year, Month, Variable, GCM)
    group_by(Year, Month, Variable, GCM) %>% 
    summarise(
      # KEY: Use the 'Mean' (daily spatial mean) column to calculate the monthly mean
      Month_Avg = mean(Mean, na.rm = TRUE), 
      N_days = n(), # Count of days in the month
      .groups = "drop"
    )
  
  return(df_monthly_result)
}


# --- Execution ---
# First, iterate over all NetCDF files using lapply
full_data_list_hist <- lapply(all_nc_files, process_nc_terra_month)

# Second, combine the list of resulting data frames into one final table
final_data_hist <- bind_rows(full_data_list_hist)

# Display the final result
print(head(final_data_hist))

###### process the Togo station data #######

library(tidyverse)
library(lubridate)
library(terra) # Required for consistency with your working environment

# =========================================================================
# === 1. GLOBAL PARAMETERS AND INITIAL READING ===
# =========================================================================

# **CONFIGURATION REQUIRED:** Set the folder containing station CSVs
PATH_TO_STATION_CSVS <- "PATH_TO_Togo_MeteoStations"

# Historical validation period
start_hist <- as.Date("1979-01-01")
end_hist <- as.Date("2014-12-31")

# --- Unified Daily Data (station_all) ---
# Find all daily CSV files (excluding quality files which start with "Datos_")
daily_files <- list.files(PATH_TO_STATION_CSVS, pattern = "\\.csv$", full.names = TRUE)
daily_files <- daily_files[!grepl("^Datos_", basename(daily_files))]

# Read all daily files and combine them
station_all <- lapply(daily_files, read_csv, show_col_types = FALSE) |>
  bind_rows() |>
  mutate(
    DATE = as.Date(DATE),
    
    # Correction: Divide by 10 to convert from tenths of mm to mm.
    PRCP = as.numeric(PRCP) / 10, 
    
    # Temperatures: Correction from 0.1°C to °C
    TMAX = as.numeric(TMAX) / 10,
    TMIN = as.numeric(TMIN) / 10,
    
    # CALCULATION OF TAS: Daily Mean Temperature
    TAS = (TMAX + TMIN) / 2
  )

# --- Unified Quality Files (quality_all) ---
# Find all quality CSV files (starting with "Datos_latlon")
quality_files <- list.files(PATH_TO_STATION_CSVS, pattern = "^Datos_latlon.*\\.csv$", full.names = TRUE)

# Read all quality files and combine them
quality_all <- lapply(quality_files, read_csv, show_col_types = FALSE) |>
  bind_rows() |>
  mutate(YEAR = as.integer(YEAR)) # Ensure Year is an integer for joining

# =========================================================================
# === 2. QUALITY AND TEMPORAL FILTERING ===
# =========================================================================

# Identify "good" years for each station based on data completeness (> 90%)
good_years_per_station <- quality_all |>
  mutate(
    # Treat NA as 0 for filtering purposes
    PRCP = coalesce(PRCP, 0), 
    TMAX = coalesce(TMAX, 0),
    TMIN = coalesce(TMIN, 0),
    # Check if either TMAX or TMIN completeness is >= 90%
    tas_ok = TMAX >= 90 | TMIN >= 90 
  ) |>
  # Keep year/station if PRCP is good OR if TAS (temp) is good
  filter(PRCP >= 90 | tas_ok) |>
  select(STATION, YEAR) # Only need to know if the year is valid

# Apply the filter to the daily data
station_valid_filtered <- station_all |>
  
  mutate(YEAR = year(DATE)) |>
  
  # Join with the table of 'good' years (inner_join discards bad years/stations)
  inner_join(good_years_per_station, by = c("STATION", "YEAR")) |>
  
  # Temporal filtering (1979-2014 period)
  filter(DATE >= start_hist, DATE <= end_hist)

# =========================================================================
# === 3. AGGREGATION TO MONTHLY AVERAGES/SUMS ===
# =========================================================================

df_monthly_validation <- station_valid_filtered |>
  
  mutate(MONTH = month(DATE, label = FALSE)) |>
  
  # Group by the combination of Station, Year, and Month
  group_by(STATION, YEAR, MONTH, LATITUDE, LONGITUDE) |>
  
  summarise(
    # PRCP (Precipitation): Use the monthly SUM
    PRCP_mm = sum(PRCP, na.rm = TRUE),
    
    # TAS, TMAX, TMIN: Use the monthly MEAN
    TAS_avg = mean(TAS, na.rm = TRUE),
    TMAX_avg = mean(TMAX, na.rm = TRUE),
    TMIN_avg = mean(TMIN, na.rm = TRUE),
    
    Days_Counted = n(), # Number of valid days counted in the month
    
    .groups = "drop" 
  ) |>
  
  # Pivot to long format to facilitate KGE calculation by variable
  pivot_longer(
    cols = c(PRCP_mm, TAS_avg, TMAX_avg, TMIN_avg),
    names_to = "Variable",
    values_to = "Obs_Monthly" # Monthly Observation Value
  )

# =========================================================================
# === 4. RENAMING AND FINAL RESULT ===
# =========================================================================

# Rename variables to match the NEX-GDDP-CMIP6 names (pr, tas, tasmax, tasmin)
df_monthly_validation <- df_monthly_validation |>
  mutate(
    Variable = case_match(
      Variable,
      "PRCP_mm" ~ "pr",
      "TAS_avg" ~ "tas",
      "TMAX_avg" ~ "tasmax",
      "TMIN_avg" ~ "tasmin"
    )
  )

print(head(df_monthly_validation))

####### KGE and data preparation #####

library(tidyverse)
if (!require("hydroGOF")) install.packages("hydroGOF") # For the KGE calculation function
library(lubridate)

# ====================================================================
# === 1. KGE FUNCTION AND INPUT DATA SETUP ===
# ====================================================================

# Helper function to calculate KGE
calculate_kge <- function(obs, sim) {
  # Clean up: Remove rows with NA in either observation or simulation
  df_clean <- na.omit(data.frame(Obs = obs, Sim = sim))
  
  if (nrow(df_clean) < 2) {
    # Cannot calculate KGE without at least 2 data pairs
    return(NA_real_) 
  }
  
  # KGE from hydroGOF returns a vector, we take the first element (the KGE value)
  kge_result <- KGE(sim = df_clean$Sim, obs = df_clean$Obs, na.rm = FALSE)[1]
  return(kge_result)
}

# Assign data frames from previous steps
# Observation data (from stations)
df_obs_monthly <- df_monthly_validation
# Downscaled data (from NetCDF)
model_data_monthly <- final_data_hist


# ====================================================================
# === 2. STATION DATA PREPARATION (REGIONAL OBSERVATION) ===
# ====================================================================

# Average the station observations to a REGIONAL level using WEIGHTING
df_obs_regional <- df_obs_monthly %>%
  
  # Calculate the weighted contribution of each station
  group_by(YEAR, MONTH, Variable) %>%
  summarise(
    # --- WEIGHTED AVERAGE CALCULATION ---
    # Weight = Days_Counted
    Obs_Regional = sum(Obs_Monthly * Days_Counted, na.rm = TRUE) / 
      sum(Days_Counted, na.rm = TRUE),
    # ------------------------------------
    .groups = "drop"
  ) %>%
  rename(Year = YEAR, Month = MONTH)


# ====================================================================
# === 3. DOWNSCALED DATA PREPARATION (REGIONAL SIMULATION) ===
# ====================================================================

sim_data_regional <- model_data_monthly %>%
  mutate(
    # For 'pr', convert from Monthly Average Daily Rate to Monthly Total
    Sim_Regional = case_when(
      Variable == "pr" ~ Month_Avg * N_days,
      # For other variables, Month_Avg is the desired monthly mean.
      TRUE ~ Month_Avg 
    )
  ) %>%
  select(Year, Month, Variable, GCM, Sim_Regional)


# ====================================================================
# === 4. MERGE AND KGE CALCULATION (Non-Bootstrap) ===
# ====================================================================

# Merge the Observation and Simulation time series
df_kge_merged <- inner_join(
  sim_data_regional,
  df_obs_regional,
  by = c("Year", "Month", "Variable")
)

# Calculate KGE for the complete time series (GCM x Variable)
df_kge_result <- df_kge_merged %>%
  
  group_by(GCM, Variable) %>%
  summarise(
    KGE = calculate_kge(obs = Obs_Regional, sim = Sim_Regional),
    N = n(), # Number of data points used
    .groups = "drop"
  )

print(df_kge_result)

# ====================================================================
# === 5. BOOTSTRAP FUNCTION PER GROUP (GCM x Variable) ===
# ====================================================================

# Define the number of Bootstrap iterations
B <- 1000 

# Function to perform Bootstrap on a time series
boot_kge_calc <- function(data, B_iterations) {
  
  kge_results <- numeric(B_iterations)
  n_samples <- nrow(data)
  
  for (i in 1:B_iterations) {
    # 1. Resampling with replacement (Bootstrap)
    boot_indices <- sample(1:n_samples, size = n_samples, replace = TRUE)
    
    # 2. Create the Bootstrap sample of the data
    boot_data <- data[boot_indices, ]
    
    # 3. Calculate KGE for this Bootstrap sample
    kge_results[i] <- calculate_kge(
      obs = boot_data$Obs_Regional, 
      sim = boot_data$Sim_Regional
    )
  }
  
  # 4. Summarize the results
  # Calculate the median and 95% confidence interval 
  kge_final <- quantile(kge_results, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  return(
    tibble(
      KGE_Median = kge_final[2],
      KGE_LCI_95 = kge_final[1], # Lower Confidence Interval (2.5%)
      KGE_UCI_95 = kge_final[3], # Upper Confidence Interval (97.5%)
      KGE_SD = sd(kge_results, na.rm = TRUE) # Standard Deviation of KGE estimates
    )
  )
}

# ====================================================================
# === 6. BOOTSTRAP APPLICATION AND FINAL RESULT ===
# ====================================================================

df_kge_bootstrap_result <- df_kge_merged %>%
  
  # Group by GCM and Variable
  group_by(GCM, Variable) %>%
  
  # Apply the Bootstrap function to each group
  do(boot_kge_calc(., B_iterations = B)) %>%
  
  ungroup()

print(df_kge_bootstrap_result)


# ====================================================================
# === 7. VISUALIZATION (KGE Plot) ===
# ====================================================================

library(ggplot2)

# --- Create the KGE Performance Plot ---

kge_plot <- df_kge_bootstrap_result %>%
  
  # Use the Median for the point position and confidence limits for the error bar
  ggplot(aes(x = GCM, y = KGE_Median, color = Variable)) +
  
  # Error Bars (95% Confidence Interval)
  geom_errorbar(
    aes(ymin = KGE_LCI_95, ymax = KGE_UCI_95),
    width = 0.2, # Width of the whisker tips
    position = position_dodge(0.7) # Separate bars by variable
  ) +
  
  # Points (KGE Median)
  geom_point(
    size = 3,
    position = position_dodge(0.7)
  ) +
  
  # Reference Line for 'Acceptable' Performance (KGE > 0.5)
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
  
  # Reference Line for 'Bad' Performance (KGE = 0.0)
  geom_hline(yintercept = 0.0, linetype = "dashed", color = "red") +
  
  # Labels and Titles
  labs(
    title = "Downscaling Performance (1979-2014) | KGE with Bootstrap",
    subtitle = "Vertical bars represent the 95% Confidence Interval",
    x = "GCM Model",
    y = "KGE (Kling-Gupta Efficiency)",
    color = "Variable"
  ) +
  
  # Y-axis Scale (adjust if your values are very different)
  scale_y_continuous(limits = c(-0.2, 1.0), breaks = seq(-0.2, 1.0, 0.2)) +
  
  # Themes and Aesthetics
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate GCM labels
    legend.position = "bottom"
  )

print(kge_plot)