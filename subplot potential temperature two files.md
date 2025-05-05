```python
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# File paths for the NetCDF files
athens_file = 'athinatwo.nc'  # File containing temperature data for Athens
thessaloniki_file = 'thessalonikitwo.nc'  # File containing temperature data for Thessaloniki

# Coordinates for Athens and Thessaloniki
athens_coords = (37.9838, 23.7275)
thessaloniki_coords = (40.6401, 22.9444)

# Constants for potential temperature calculation
Rd = 287.05  # Specific gas constant for dry air (J/kg·K)
cp = 1005  # Specific heat capacity at constant pressure (J/kg·K)
p0 = 1000  # Reference pressure level (hPa)

def potential_temperature(T, P):
    return T * (p0 / P) ** (Rd / cp)

def process_file(file, location_coords):
    ds = nc.Dataset(file, 'r')
    
    latitudes = ds.variables['latitude'][:]
    longitudes = ds.variables['longitude'][:]
    pressure_levels = ds.variables['pressure_level'][:]
    time_var = ds.variables['valid_time'][:]
    
    time_units = ds.variables['valid_time'].units
    time_calendar = getattr(ds.variables['valid_time'], 'calendar', 'standard')
    time_values = nc.num2date(time_var, units=time_units, calendar=time_calendar)
    
    def find_nearest_idx(array, value):
        return np.abs(array - value).argmin()
    
    lat_idx = find_nearest_idx(latitudes, location_coords[0])
    lon_idx = find_nearest_idx(longitudes, location_coords[1])
    
    temperature = ds.variables['t'][:, :, lat_idx, lon_idx] - 273.15  # Convert K to Celsius
    geopotential_height = (ds.variables['z'][:, :, lat_idx, lon_idx] / 9.81) / 1000  # Convert to km
    
    # Compute potential temperature
    potential_temp = potential_temperature(temperature + 273.15, pressure_levels[None, :]) - 273.15  # Convert back to Celsius
    
    ds.close()
    
    return potential_temp, geopotential_height, time_values, pressure_levels

athens_potential_temp, athens_geopotential, athens_time_values, pressure_levels = process_file(athens_file, athens_coords)
thessaloniki_potential_temp, thessaloniki_geopotential, thessaloniki_time_values, _ = process_file(thessaloniki_file, thessaloniki_coords)

def filter_season_indices(time_values, months):
    return [i for i, t in enumerate(time_values) if t.month in months]

winter_months = [12, 1, 2]
spring_months = [3, 4, 5]
summer_months = [6, 7, 8]
autumn_months = [9, 10, 11]

athens_winter_indices = filter_season_indices(athens_time_values, winter_months)
thessaloniki_winter_indices = filter_season_indices(thessaloniki_time_values, winter_months)
athens_spring_indices = filter_season_indices(athens_time_values, spring_months)
thessaloniki_spring_indices = filter_season_indices(thessaloniki_time_values, spring_months)
athens_summer_indices = filter_season_indices(athens_time_values, summer_months)
thessaloniki_summer_indices = filter_season_indices(thessaloniki_time_values, summer_months)
athens_autumn_indices = filter_season_indices(athens_time_values, autumn_months)
thessaloniki_autumn_indices = filter_season_indices(thessaloniki_time_values, autumn_months)

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

def calculate_average_profile(data, indices):
    return np.mean(data[indices, :], axis=0) if indices else None

fig, axs = plt.subplots(2, 2, figsize=(15, 12))

def plot_season_profiles(ax, season, athens_time_indices, thessaloniki_time_indices, athens_data, athens_geopotential, thessaloniki_data, thessaloniki_geopotential):
    for time_label, indices in athens_time_indices.items():
        avg_profile = calculate_average_profile(athens_data, indices)
        if avg_profile is not None:
            ax.plot(avg_profile, athens_geopotential[0, :], label=f'Athens - {time_label}', linestyle='--')
    
    for time_label, indices in thessaloniki_time_indices.items():
        avg_profile = calculate_average_profile(thessaloniki_data, indices)
        if avg_profile is not None:
            ax.plot(avg_profile, thessaloniki_geopotential[0, :], label=f'Thessaloniki - {time_label}', linestyle='-')
    
    ax.set_ylabel('Geopotential Height (km)')
    ax.set_xlabel('Potential Temperature (°C)')
    ax.set_title(f'Average Vertical Profile of Potential Temperature in {season}')
    ax.legend()
    ax.set_xlim(-5, 40)
    ax.set_ylim(0, 3)

plot_season_profiles(axs[0, 0], 'Winter', athens_winter_time_indices, thessaloniki_winter_time_indices, athens_potential_temp, athens_geopotential, thessaloniki_potential_temp, thessaloniki_geopotential)
plot_season_profiles(axs[0, 1], 'Spring', athens_spring_time_indices, thessaloniki_spring_time_indices, athens_potential_temp, athens_geopotential, thessaloniki_potential_temp, thessaloniki_geopotential)
plot_season_profiles(axs[1, 0], 'Summer', athens_summer_time_indices, thessaloniki_summer_time_indices, athens_potential_temp, athens_geopotential, thessaloniki_potential_temp, thessaloniki_geopotential)
plot_season_profiles(axs[1, 1], 'Autumn', athens_autumn_time_indices, thessaloniki_autumn_time_indices, athens_potential_temp, athens_geopotential, thessaloniki_potential_temp, thessaloniki_geopotential)

plt.tight_layout()
plt.show()
```


    
![png](output_0_0.png)
    



```python

```
