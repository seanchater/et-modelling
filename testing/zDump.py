from cmath import inf
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
import outputs as out
import PARAMS as parm
import os

vc = ""
ndvi = ""


# **effective_leaf_area_index**************************************************
# constants or predefined:
""" nd_min = np.nanmin(ndvi)
nd_max = np.nanmax(ndvi)
vc_pow = 0.7
vc_min = np.nanmin(vc)
vc_max = np.nanmax(vc)
lai_pow = -0.45 """



""" def exception(test):
    try:
        float(test)
        return True
    except:
        return False

lai_pow = 0.0
vc_pow = 0.0

def setPowLAI(value):
    global lai_pow
    lai_pow = value

def getPowLAI():
    global lai_pow
    return lai_pow

in_ = input("LAI pow: ")
setPowLAI([lambda:-0.45, lambda: float(in_)][exception(in_)]())
print(getPowLAI())
print(exception(in_)) """


def clipRast(ras):
    '''
    (GDAL raster, GDAL extent)
    '''
    try:
    # newRas = gdal.Warp("gdalwarp -cutline G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/3318DD_extent/3318DD.shp -crop_to_cutline G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/DEM_Derviatives/aspect.tif cropped_output.tif")
    # newRas = gdal.Warp(gdal.WarpOptions(" -cutline G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/3318DD_extent/3318DD.shp -crop_to_cutline G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/DEM_Derviatives/aspect.tif cropped_output.tif"))
        gdal.Warp(os.path.join(), ras, cutlineDSName="G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/3318DD_extent/3318DD.shp", cropToCutline=ras)
        return True
    except:
        return False
    
# def clipFail():
    #IN ORDER TO CLIP BY EXTENT EVERY IMAGE
    """ gt1=ras.GetGeoTransform()
    gt2=extent.GetGeoTransform()
    if gt1[0] < gt2[0]: #CONDITIONAL TO SELECT THE CORRECT ORIGIN
        gt3=gt2[0]
    else:
        gt3=gt1[0]
    if gt1[3] < gt2[3]:
        gt4=gt1[3]
    else:
        gt4=gt2[3]
    xOrigin = gt3
    yOrigin = gt4
    pixelWidth = gt1[1]
    pixelHeight = gt1[5]

    r1 = [gt1[0], gt1[3],gt1[0] + (gt1[1] * ras.RasterXSize), gt1[3] + (gt1[5] * ras.RasterYSize)]
    r2 = [gt2[0], gt2[3],gt2[0] + (gt2[1] * extent.RasterXSize), gt2[3] + (gt2[5] * extent.RasterYSize)]
    intersection = [max(r1[0], r2[0]), min(r1[1], r2[1]), min(r1[2], r2[2]), max(r1[3], r2[3])]

    xmin = intersection[0]
    xmax = intersection[2]
    ymin = intersection[3]
    ymax = intersection[1]

    print("x {min, max} : ", xmin, xmax, "\ny {min, max} : ", ymin, ymax)

    # Specify offset and rows and columns to read
    xoff = int((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1
    srs=ras.GetProjectionRef() #necessary to export with SRS

    ras = gdal.Translate("G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/tst/cliptst.tif", ras, srcWin = [xoff, yoff, xcount, ycount])
    ras.ReadAsArray()
    ras = np.round(ras, 4)
    return ras """
    """ target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((xmin, pixelWidth, 0,ymax, 0, pixelHeight,))
    img1 = ras.GetRasterBand(1).ReadAsArray(xoff, yoff, xcount, ycount)  
    img2 = extent.GetRasterBand(1).ReadAsArray(xoff, yoff, xcount, ycount)

    return target_ds """

""" def clipRast(ras):
    '''
    (GDAL raster, GDAL extent)
    '''
    try:
        gdal.Warp((input_file_path), ras, cutlineDSName="G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/3318DD_extent/3318DD.shp", cropToCutline=ras)
        return True
    except:
        return False """

# gdalwarp -overwrite -to SRC_METHOD=NO_GEOTRANSFORM -to DST_METHOD=NO_GEOTRANSFORM -te 220 60 1160 690 -cutline cutline.csv in.png out.tif##########################

dest_ext = gdal.Open("G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/tst/extent.tif")
ext = dest_ext.GetRasterBand(1).ReadAsArray()
# ext = np.round(ext, 4)
# ext[ext == -inf] = np.nan
dest_slope = gdal.Open("G:/My Drive/Stellenbosch/2022/716/ET/modeling/Data/in_/DEM_Derviatives/aspect.tif")
aspect_np = dest_slope.GetRasterBand(1).ReadAsArray()
aspect_np = np.round(aspect_np, 4)
aspect_np[aspect_np == inf] = np.nan
aspect_deg = clipRast(dest_slope)


print(ext)
print("aspect_np ", aspect_np, "\nmin: ", np.nanmin(aspect_np), "\nmax: ", np.nanmax(aspect_np))
print("aspect_deg ", aspect_deg, "\nmin: ", np.nanmin(aspect_deg), "\nmax: ", np.nanmax(aspect_deg))


############################################################################################## OG INPUTS
""" 

dest_lst = gdal.Open(par.getFilePathIN("lst"))
lst = dest_lst.GetRasterBand(1).ReadAsArray()
lst = np.round(lst, 4)
lst[lst == -inf] = np.nan
print("LST ", lst, "\nmin: ", np.nanmin(lst), "\nmax: ", np.nanmax(lst)) #tst

dest_albedo = gdal.Open(par.getFilePathIN("albedo"))
r0 = dest_albedo.GetRasterBand(1).ReadAsArray()
r0 = np.round(r0, 4)
r0[r0 == inf] = np.nan
print("R0 ", r0, "\nmin: ", np.nanmin(r0), "\nmax: ", np.nanmax(r0)) #tst

dest_ndvi = gdal.Open(par.getFilePathIN("ndvi"))
ndvi = dest_ndvi.GetRasterBand(1).ReadAsArray()
ndvi = np.round(ndvi, 4)
ndvi[ndvi == -inf] = np.nan
print("NDVI ", ndvi, "\nmin: ", np.nanmin(ndvi), "\nmax: ", np.nanmax(ndvi)) #tst

desttime = gdal.Open(par.getFilePathIN("time"))
dtime = desttime.GetRasterBand(1).ReadAsArray()
dtime = np.round(dtime, 4)
dtime[dtime == -inf] = np.nan
print("dtime ", dtime, "\nmin: ", np.nanmin(dtime), "\nmax: ", np.nanmax(dtime)) #tst

dest_lat = gdal.Open(par.getFilePathIN("lat"))
lat_deg = dest_lat.GetRasterBand(1).ReadAsArray()
print("Lat ", lat_deg, "\nmin: ", np.nanmin(lat_deg), "\nmax: ", np.nanmax(lat_deg)) #tst

dest_lon = gdal.Open(par.getFilePathIN("lon"))
lon_deg = dest_lon.GetRasterBand(1).ReadAsArray()
print("Lon ", lon_deg, "\nmin: ", np.nanmin(lon_deg), "\nmax: ", np.nanmax(lon_deg)) #tst

dest_dem = gdal.Open(par.getFilePathIN("dem"))
z = dest_dem.GetRasterBand(1).ReadAsArray()
#z[np.isnan(lst)] = np.nan
print("Z ", z, "\nmin: ", np.nanmin(z), "\nmax: ", np.nanmax(z)) #tst

dest_slope = gdal.Open(par.getFilePathIN("slope"))
slope_deg = dest_slope.GetRasterBand(1).ReadAsArray()
# slope_deg = np.round(slope_deg, 5)
# slope_deg[slope_deg == inf] = np.nan
print("Slope ", slope_deg, "\nmin: ", np.nanmin(slope_deg), "\nmax: ", np.nanmax(slope_deg)) #tst

dest_aspect = gdal.Open(par.getFilePathIN("aspect"))
aspect_deg = dest_aspect.GetRasterBand(1).ReadAsArray()
# aspect_deg = np.round(aspect_deg, 8)
# aspect_deg[aspect_deg == inf] = np.nan
print("Aspect ", aspect_deg, "\nmin: ", np.nanmin(aspect_deg), "\nmax: ", np.nanmax(aspect_deg)) #tst

dest_lm = gdal.Open(par.getFilePathIN("landMask"))
land_mask = dest_lm.GetRasterBand(1).ReadAsArray()
#land_mask[np.isnan(lst)] = np.nan
print("LandMask ", land_mask, "\nmin: ", np.nanmin(land_mask), "\nmax: ", np.nanmax(land_mask)) #tst

#dest_bulk = gdal.Open(par.getFilePathIN("bulk"))
#bulk = dest_bulk.GetRasterBand(1).ReadAsArray()

dest_maxobs = gdal.Open(par.getFilePathIN("maxObs"))
z_obst_max = dest_maxobs.GetRasterBand(1).ReadAsArray()
#z_obst_max[np.isnan(lst)] = np.nan
print("ZobsMax ", z_obst_max, "\nmin: ", np.nanmin(z_obst_max), "\nmax: ", np.nanmax(z_obst_max)) #tst

dest_pairsea24 = gdal.Open(par.getFilePathIN("pair_24_0"))
p_air_0_24 = dest_pairsea24.GetRasterBand(1).ReadAsArray()
p_air_0_24 = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_24)
p_air_0_24[np.isnan(lst)] = np.nan

dest_pairseainst = gdal.Open(par.getFilePathIN("pair_inst_0"))
p_air_0_i = dest_pairseainst.GetRasterBand(1).ReadAsArray()
p_air_0_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_i)
p_air_0_i[np.isnan(lst)] = np.nan

dest_pairinst = gdal.Open(par.getFilePathIN("pair_inst"))
p_air_i = dest_pairinst.GetRasterBand(1).ReadAsArray()
p_air_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_i)
p_air_i[np.isnan(lst)] = np.nan

dest_precip = gdal.Open(par.getFilePathIN("pre"))
P_24 = dest_precip.GetRasterBand(1).ReadAsArray()
P_24[np.isnan(lst)] = np.nan

dest_hum24 = gdal.Open(par.getFilePathIN("hum_24"))
qv_24 = dest_hum24.GetRasterBand(1).ReadAsArray()
qv_24[np.isnan(lst)] = np.nan

dest_huminst = gdal.Open(par.getFilePathIN("hum_inst"))
qv_i = dest_huminst.GetRasterBand(1).ReadAsArray()
qv_i[np.isnan(lst)] = np.nan

dest_tair24 = gdal.Open(par.getFilePathIN("tair_24"))
t_air_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
#t_air_24 = ETLook.meteo.disaggregate_air_temperature_daily(t_air_24_coarse, z, z_coarse, lapse)
t_air_24[np.isnan(lst)] = np.nan

dest_tair24 = gdal.Open(par.getFilePathIN("tair_max_24"))
t_air_max_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
t_air_max_24[np.isnan(lst)] = np.nan

dest_tair24 = gdal.Open(par.getFilePathIN("tair_min_24"))
t_air_min_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
t_air_min_24[np.isnan(lst)] = np.nan

dest_tairinst = gdal.Open(par.getFilePathIN("tair_inst"))
t_air_i = dest_tairinst.GetRasterBand(1).ReadAsArray()
t_air_i[np.isnan(lst)] = np.nan

dest_tairamp = gdal.Open(par.getFilePathIN("tair_amp"))
t_amp_year = dest_tairamp.GetRasterBand(1).ReadAsArray()
t_amp_year[np.isnan(lst)] = np.nan

dest_wind24 = gdal.Open(par.getFilePathIN("wind_24"))
u_24 = dest_wind24.GetRasterBand(1).ReadAsArray()
u_24[np.isnan(lst)] = np.nan

dest_windinst = gdal.Open(par.getFilePathIN("wind_inst"))
u_i = dest_windinst.GetRasterBand(1).ReadAsArray()
u_i[np.isnan(lst)] = np.nan

dest_watcol = gdal.Open(par.getFilePathIN("watCol_inst"))
wv_i = dest_watcol.GetRasterBand(1).ReadAsArray()
wv_i[np.isnan(lst)] = np.nan

dest_trans = gdal.Open(par.getFilePathIN("trans_24"))
trans_24 = dest_trans.GetRasterBand(1).ReadAsArray()
#trans_24[np.isnan(lst)] = np.nan
print("Trans24 ", trans_24, "\nmin: ", np.nanmin(trans_24), "\nmax: ", np.nanmax(trans_24)) #tst

 """