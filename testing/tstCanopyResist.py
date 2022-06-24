import os
from weakref import ref
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
import ref_solar_radiation as sr
import ref_clear_sky_radiation as csr
import datetime as dt

input_folder = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data"
output_folder = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data\output"

print("Starting with:\nInput Folder  - ", input_folder, "\nOutput Folder - ", output_folder)

# Input folder Date
input_folder_date = os.path.join(input_folder)
# Output folder Date
output_folder_date = os.path.join(output_folder)
if not os.path.exists(output_folder_date):
    os.makedirs(output_folder_date)

date = dt.datetime(2021, 1, 1)
doy = int(date.strftime("%j"))

# lat
lat_filename = os.path.join(input_folder_date, r"latlon\lat_v2.tif")
dest_lat = gdal.Open(lat_filename)
print("\n\n", dest_lat)
print("X ", dest_lat.RasterXSize)
print("Y ", dest_lat.RasterYSize)
print(type(dest_lat))
lat_deg = dest_lat.ReadAsArray()
# lat_deg[np.isnan(lst)] = np.nan

# slope 
slope_filename = os.path.join(input_folder_date, r"DEM_derivatives_2\Slope_AOI.tif")
dest_slope = gdal.Open(slope_filename)
print("\n\n", dest_slope)
print("X ", dest_slope.RasterXSize)
print("Y ", dest_slope.RasterYSize)
print(type(dest_slope))
slope_deg = dest_slope.GetRasterBand(1).ReadAsArray()
# slope_deg[np.isnan(lst)] = np.nan

aspect_filename = os.path.join(input_folder_date, r"DEM_derivatives_2\Aspect_AOI.tif")
dest_aspect = gdal.Open(aspect_filename)
print("\n\n", dest_aspect)
print("X ", dest_aspect.RasterXSize)
print("Y ", dest_aspect.RasterYSize)
print(type(dest_aspect))
aspect_deg = dest_aspect.ReadAsArray()
# aspect_deg[np.isnan(lst)] = np.nan

# get ref system for saving later
# geo_ex = lat_deg.GetGeoTransform()
# proj_ex = lat_deg.GetProjection()

##### TRANSPIRATION COMPONENT #################################
# constants or predefined:
diffusion_slope = -1.33
diffusion_intercept = 1.15
t_opt = 25 # optimal temperature for plant growth
t_min = 0 # minimal temperature for plant growth
t_max = 50 # maximal temperature for plant growth
vpd_slope = -0.3
rs_min = 70
rcan_max = 1000000

# ** atmospheric canopy resistance ********************************
iesd = sr.inverse_earth_sun_distance(doy)
sc = sr.seasonal_correction(doy)
day_angle = csr.day_angle(doy)
decl = csr.declination(doy)

lat = sr.latitude_rad(lat_deg)
slope = sr.slope_rad(slope_deg)
aspect = sr.aspect_rad(aspect_deg)
ra_24_toa = sr.daily_solar_radiation_toa(sc, decl, iesd, lat, slope, aspect)
ws = sr.sunset_hour_angle(lat, decl)
ra_24_toa_flat = sr.daily_solar_radiation_toa_flat(decl, iesd, lat, ws)
# diffusion_index = sr.diffusion_index(trans_24, diffusion_slope, diffusion_intercept)
