import cdsapi
import os

c = cdsapi.Client()
output_folder = '/home/david/Desktop/ERA5_Pmm'
os.makedirs(output_folder, exist_ok=True)

# Variables for precipitation
variables_precip = ['total_precipitation']

years_list = [str(year) for year in range(1940, 2025)]

# Download precipitation data - one file per month with all years
for month in range(1, 13):
    month_str = f"{month:02d}"
    output_file = os.path.join(output_folder,
f'era5_precip_month_{month_str}.nc')

    c.retrieve(
        'reanalysis-era5-land-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': variables_precip,
            'year': years_list,  # All years in one request
            'month': month_str,  # Single month
            'time': '00:00',
            'data_format': 'netcdf',
            # 'area': [90, -180, -90, 180],  # Optional: explicitly
set global bounds
        },
        output_file)
    print(f"Completed download for month {month_str}")

print("All downloads completed!")
