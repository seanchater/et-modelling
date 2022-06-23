import os
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
import leaf
import Processing_Functions as PF
#import ETLook_
#import Functions_ as PF
import outputs as out

input_folder = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data\Indecies\Indecies"
output_folder = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data\output"

print("Starting with:\nInput Folder  - ", input_folder, "\nOutput Folder - ", output_folder)

# Input folder Date
input_folder_date = os.path.join(input_folder)
# Output folder Date
output_folder_date = os.path.join(output_folder)
if not os.path.exists(output_folder_date):
    os.makedirs(output_folder_date)

NDVI_filename = os.path.join(input_folder_date, "Indeces_3318DD_Jan.tif")
vc_filename = os.path.join(output_folder_date, "vc_.tif")
lai_filename = os.path.join(output_folder_date, "LAI_.tif")
lai_eff_filename= os.path.join(output_folder_date, "LAI_eff_.tif")


dest_ndvi = gdal.Open(NDVI_filename)

b1ndvi = dest_ndvi.GetRasterBand(1)
d1ndvi = b1ndvi.ReadAsArray()

plt.imshow(d1ndvi)
plt.show()

ndvi = dest_ndvi.GetRasterBand(1).ReadAsArray()

plt.hist(ndvi, bins=5)
plt.show()

print(ndvi)
ndvi[ndvi == -9999] = np.nan
print(ndvi)
#ndvi[np.isnan(lst)] = np.nan

# example file
geo_ex = dest_ndvi.GetGeoTransform()
proj_ex = dest_ndvi.GetProjection()

# Create QC array
QC = np.ones(ndvi.shape)
QC[np.isnan(ndvi)] = np.nan

# page 31 flow diagram

# **effective_leaf_area_index**************************************************
# constants or predefined:
nd_min = 0.125
nd_max = 0.8
vc_pow = 0.7
vc_min = 0
vc_max = 0.9677324224821418
lai_pow = -0.45


vc = leaf.vegetation_cover(ndvi, nd_min, nd_max, vc_pow)
lai = leaf.leaf_area_index(vc, vc_min, vc_max, lai_pow)
lai_eff = leaf.effective_leaf_area_index(lai)

print("\n\n", vc)

vc[np.isnan(QC)] = np.nan
lai[np.isnan(QC)] = np.nan
lai_eff[np.isnan(QC)] = np.nan
if vc.any():
    PF.Save_as_tiff(vc_filename, vc, geo_ex, proj_ex)
    print("vc Saved...", vc_filename)
if lai.any():
    PF.Save_as_tiff(lai_filename, lai, geo_ex, proj_ex)
    print("lai Saved...", lai_filename)
if lai_eff.any():
    PF.Save_as_tiff(lai_eff_filename, lai_eff, geo_ex, proj_ex)
    print("lai_eff Saved...", lai_eff_filename)

