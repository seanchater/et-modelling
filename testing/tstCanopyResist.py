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

lat_filename = os.path.join(input_folder_date, r"\latlon\lat.tif")
lat_deg = gdal.Open(lat_filename)
lat_deg[np.isnan(lst)] = np.nan

# get ref system for saving later
geo_ex = lat_deg.GetGeoTransform()
proj_ex = lat_deg.GetProjection()

##### TRANSPIRATION COMPONENT #################################

# ** atmospheric canopy resistance ********************************
iesd = sr.inverse_earth_sun_distance(doy)
sc = sr.seasonal_correction(doy)
day_angle = csr.day_angle(doy)
decl = csr.declination(doy)

lat = sr.latitude_rad(lat_deg)
slope = sr.slope_rad(slope_deg)
aspect = sr.aspect_rad(aspect_deg)

print(iesd)