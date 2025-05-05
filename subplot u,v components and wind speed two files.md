```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# File paths for the NetCDF files
athens_file = 'athinatwo.nc'  # File containing wind data for Athens
thessaloniki_file = 'thessalonikitwo.nc'  # File containing wind data for Thessaloniki

# Coordinates for Athens and Thessaloniki
athens_coords = (37.9838, 23.7275)
thessaloniki_coords = (40.6401, 22.9444)

def process_file(file, location_coords):
    # Open the NetCDF file
    ds = nc.Dataset(file, 'r')

    # Extract latitude, longitude, pressure levels, and time arrays
    latitudes = ds.variables['latitude'][:]  # Assuming 'latitude' is the name of the lat variable
    longitudes = ds.variables['longitude'][:]  # Assuming 'longitude' is the name of the lon variable
    pressure_levels = ds.variables['pressure_level'][:]  # Assuming 'pressure_level' is the name
    time_var = ds.variables['valid_time'][:]  # Assuming 'valid_time' is the name of the time variable

    # Convert time units
    time_units = ds.variables['valid_time'].units
    time_calendar = getattr(ds.variables['valid_time'], 'calendar', 'standard')
    time_values = nc.num2date(time_var, units=time_units, calendar=time_calendar)

    # Find the nearest grid point for the location
    def find_nearest_idx(array, value):
        return np.abs(array - value).argmin()

    lat_idx = find_nearest_idx(latitudes, location_coords[0])
    lon_idx = find_nearest_idx(longitudes, location_coords[1])

    # Extract the u and v wind components
    u_component = ds.variables['u'][:]  # Shape: (time, pressure_levels, lat, lon)
    v_component = ds.variables['v'][:]  # Shape: (time, pressure_levels, lat, lon)

    # Calculate wind speed (sqrt(u^2 + v^2))
    wind_speed = np.sqrt(u_component**2 + v_component**2)

    # Extract the geopotential height or related variable (z/g)
    geopotential_height = ds.variables['z'][:]  # Assuming 'z' is the name of the variable

    # Extract u, v, wind speed and geopotential height time series for the location
    u_profile = u_component[:, :, lat_idx, lon_idx]  # Shape: (time, pressure_levels)
    v_profile = v_component[:, :, lat_idx, lon_idx]  # Shape: (time, pressure_levels)
    wind_speed_profile = wind_speed[:, :, lat_idx, lon_idx]  # Shape: (time, pressure_levels)
    geopotential_profile = (geopotential_height[:, :, lat_idx, lon_idx] / 9.81) / 1000  # Convert z to z/g (km)

    # Close the NetCDF file
    ds.close()

    return u_profile, v_profile, wind_speed_profile, geopotential_profile, time_values, pressure_levels

# Process files for Athens and Thessaloniki
athens_u, athens_v, athens_wind_speed, athens_geopotential, athens_time_values, pressure_levels = process_file(athens_file, athens_coords)
thessaloniki_u, thessaloniki_v, thessaloniki_wind_speed, thessaloniki_geopotential, thessaloniki_time_values, _ = process_file(thessaloniki_file, thessaloniki_coords)

# Filter time values for each season
def filter_season_indices(time_values, months):
    return [i for i, t in enumerate(time_values) if t.month in months]

# Define the months for each season
winter_months = [12, 1, 2]
spring_months = [3, 4, 5]
summer_months = [6, 7, 8]
autumn_months = [9, 10, 11]

# Get indices for each season
athens_winter_indices = filter_season_indices(athens_time_values, winter_months)
thessaloniki_winter_indices = filter_season_indices(thessaloniki_time_values, winter_months)

athens_spring_indices = filter_season_indices(athens_time_values, spring_months)
thessaloniki_spring_indices = filter_season_indices(thessaloniki_time_values, spring_months)

athens_summer_indices = filter_season_indices(athens_time_values, summer_months)
thessaloniki_summer_indices = filter_season_indices(thessaloniki_time_values, summer_months)

athens_autumn_indices = filter_season_indices(athens_time_values, autumn_months)
thessaloniki_autumn_indices = filter_season_indices(thessaloniki_time_values, autumn_months)

# Extract `00:00` and `12:00` indices for each season
def extract_time_indices(time_values, season_indices):
    return {
        '00:00': [i for i in season_indices if time_values[i].hour == 0],
        '12:00': [i for i in season_indices if time_values[i].hour == 12]
    }

athens_winter_time_indices = extract_time_indices(athens_time_values, athens_winter_indices)
thessaloniki_winter_time_indices = extract_time_indices(thessaloniki_time_values, thessaloniki_winter_indices)

athens_spring_time_indices = extract_time_indices(athens_time_values, athens_spring_indices)
thessaloniki_spring_time_indices = extract_time_indices(thessaloniki_time_values, thessaloniki_spring_indices)

athens_summer_time_indices = extract_time_indices(athens_time_values, athens_summer_indices)
thessaloniki_summer_time_indices = extract_time_indices(thessaloniki_time_values, thessaloniki_summer_indices)

athens_autumn_time_indices = extract_time_indices(athens_time_values, athens_autumn_indices)
thessaloniki_autumn_time_indices = extract_time_indices(thessaloniki_time_values, thessaloniki_autumn_indices)

# Function to calculate the average vertical profile for u, v, or wind speed
def calculate_average_profile(data, indices):
    return np.mean(data[indices, :], axis=0) if indices else None

# Create subplots for each season
fig, axs = plt.subplots(4, 3, figsize=(18, 20))

def plot_season_profiles(ax, season, athens_time_indices, thessaloniki_time_indices, athens_data, thessaloniki_data, data_label, y_label, title_label):
    for time_label, indices in athens_time_indices.items():
        avg_profile = calculate_average_profile(athens_data, indices)
        if avg_profile is not None:
            ax.plot(avg_profile, athens_geopotential[0, :], label=f'Athens - {time_label}', linestyle='--')

    for time_label, indices in thessaloniki_time_indices.items():
        avg_profile = calculate_average_profile(thessaloniki_data, indices)
        if avg_profile is not None:
            ax.plot(avg_profile, thessaloniki_geopotential[0, :], label=f'Thessaloniki - {time_label}', linestyle='-')

    ax.set_ylabel('Geopotential Height (km)')
    ax.set_xlabel(data_label)
    ax.set_title(f'Average Vertical Profile of {title_label} in {season}')
    ax.legend()
    ax.set_xlim(-30, 30)  # Adjust x-axis for u and v components, and for wind speed
    ax.set_ylim(0, 3)  # Zoom in on y-axis to show height range (0 to 3 km)

# Plot for Winter - u component
plot_season_profiles(axs[0, 0], 'Winter', athens_winter_time_indices, thessaloniki_winter_time_indices, athens_u, thessaloniki_u, 'u component of Wind (m/s)', 'u component (m/s)', 'u component of Wind')

# Plot for Winter - v component
plot_season_profiles(axs[0, 1], 'Winter', athens_winter_time_indices, thessaloniki_winter_time_indices, athens_v, thessaloniki_v, 'v component of Wind (m/s)', 'v component (m/s)', 'v component of Wind')

# Plot for Winter - Wind speed
plot_season_profiles(axs[0, 2], 'Winter', athens_winter_time_indices, thessaloniki_winter_time_indices, athens_wind_speed, thessaloniki_wind_speed, 'Wind Speed (m/s)', 'Wind Speed (m/s)', 'Wind Speed')

# Plot for Spring - u component
plot_season_profiles(axs[1, 0], 'Spring', athens_spring_time_indices, thessaloniki_spring_time_indices, athens_u, thessaloniki_u, 'u component of Wind (m/s)', 'u component (m/s)', 'u component of Wind')

# Plot for Spring - v component
plot_season_profiles(axs[1, 1], 'Spring', athens_spring_time_indices, thessaloniki_spring_time_indices, athens_v, thessaloniki_v, 'v component of Wind (m/s)', 'v component (m/s)', 'v component of Wind')

# Plot for Spring - Wind speed
plot_season_profiles(axs[1, 2], 'Spring', athens_spring_time_indices, thessaloniki_spring_time_indices, athens_wind_speed, thessaloniki_wind_speed, 'Wind Speed (m/s)', 'Wind Speed (m/s)', 'Wind Speed')

# Plot for Summer - u component
plot_season_profiles(axs[2, 0], 'Summer', athens_summer_time_indices, thessaloniki_summer_time_indices, athens_u, thessaloniki_u, 'u component of Wind (m/s)', 'u component (m/s)', 'u component of Wind')

# Plot for Summer - v component
plot_season_profiles(axs[2, 1], 'Summer', athens_summer_time_indices, thessaloniki_summer_time_indices, athens_v, thessaloniki_v, 'v component of Wind (m/s)', 'v component (m/s)', 'v component of Wind')

# Plot for Summer - Wind speed
plot_season_profiles(axs[2, 2], 'Summer', athens_summer_time_indices, thessaloniki_summer_time_indices, athens_wind_speed, thessaloniki_wind_speed, 'Wind Speed (m/s)', 'Wind Speed (m/s)', 'Wind Speed')

# Plot for Autumn - u component
plot_season_profiles(axs[3, 0], 'Autumn', athens_autumn_time_indices, thessaloniki_autumn_time_indices, athens_u, thessaloniki_u, 'u component of Wind (m/s)', 'u component (m/s)', 'u component of Wind')

# Plot for Autumn - v component
plot_season_profiles(axs[3, 1], 'Autumn', athens_autumn_time_indices, thessaloniki_autumn_time_indices, athens_v, thessaloniki_v, 'v component of Wind (m/s)', 'v component (m/s)', 'v component of Wind')

# Plot for Autumn - Wind speed
plot_season_profiles(axs[3, 2], 'Autumn', athens_autumn_time_indices, thessaloniki_autumn_time_indices, athens_wind_speed, thessaloniki_wind_speed, 'Wind Speed (m/s)', 'Wind Speed (m/s)', 'Wind Speed')

# Adjust layout and show the plot
plt.tight_layout()
plt.show()

```


    
![png](output_0_0.png)
    



```python

```
