from cmath import inf
import os
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
import leaf
import Processing_Functions as PF
#import ETLook_
#import Functions_ as PF
import outputs as out

input_folder = "G:\\My Drive\\Stellenbosch\\2022\\716\\ET\\modeling\\Data"
output_folder = "G:\\My Drive\\Stellenbosch\\2022\\716\\ET\\modeling\\Data\\out_"

print("Starting with:\nInput Folder  - ", input_folder, "\nOutput Folder - ", output_folder)

# Input folder Date
input_folder_date = os.path.join(input_folder)
# Output folder Date
output_folder_date = os.path.join(output_folder)
if not os.path.exists(output_folder_date):
    os.makedirs(output_folder_date)

NDVI_filename = os.path.join(input_folder_date, "Indecies\\Indeces_3318DD_Jan.tif")

vc_filename = os.path.join(output_folder_date, "vc_2.tif")
lai_filename = os.path.join(output_folder_date, "LAI_2.tif")
lai_eff_filename = os.path.join(output_folder_date, "LAI_eff_2.tif")

dest_ndvi = gdal.Open(NDVI_filename)
ndvi = dest_ndvi.GetRasterBand(1).ReadAsArray()
ndvi = ndvi.round(4)
print("NDVI: ", np.size(ndvi))
ndvi[ndvi == -inf] = np.nan


# example file
geo_ex = dest_ndvi.GetGeoTransform()
proj_ex = dest_ndvi.GetProjection()

# Create QC array
QC = np.ones(ndvi.shape)
#QC[np.isnan(lst)] = np.nan

# page 31 flow diagram

# **effective_leaf_area_index**************************************************
# constants or predefined:

vc_pow = 0.7
lai_pow = -0.45

nd_min = np.nanmin(ndvi)#0.125
nd_max = np.nanmax(ndvi)#0.8

vc = leaf.vegetation_cover(ndvi, nd_min, nd_max, vc_pow)
vc_min = np.nanmin(vc)#0
vc_max = np.nanmax(vc)#0.9677324224821418

lai = leaf.leaf_area_index(vc, vc_min, vc_max, lai_pow)
lai_eff = leaf.effective_leaf_area_index(lai)


plt.imshow(vc)
plt.show()
plt.imshow(lai_eff)
plt.show()
plt.hist(lai_eff, bins=50)
plt.show()

print("\n\n", ndvi)
print("\n\n", vc)
print("\n\n", lai_eff)

#vc[np.isnan(QC)] = np.nan
#lai[np.isnan(QC)] = np.nan
#lai_eff[np.isnan(QC)] = np.nan
# Saving files
if vc.any():
    PF.Save_as_tiff(vc_filename, vc, geo_ex, proj_ex)
    print("vc Saved...", vc_filename)
if lai.any():
    PF.Save_as_tiff(lai_filename, lai, geo_ex, proj_ex)
    print("lai Saved...", lai_filename)
if lai_eff.any():
    PF.Save_as_tiff(lai_eff_filename, lai_eff, geo_ex, proj_ex)
    print("lai_eff Saved...", lai_eff_filename)

    print(vc.dtype)
    print(lai.dtype)
    print(lai_eff.dtype)