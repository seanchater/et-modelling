from cmath import inf
import os
import sys
import datetime
# from osgeo import osr
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np
import outputs as out
import PARAMS as parm
# import PARAMclip as cparm
import tkinter as tk
from tkinter import filedialog
import csv
import math
import zeditETLook, solar_radiation, clear_sky_radiation, meteo, radiation, evapotranspiration, soil_moisture, leaf, stress, resistance, roughness, neutral, unstable, outputs
# from Functions import Processing_Functions as PF

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
print(parent)
# import ETLook as ET
import Functions as PF

#np.set_printoptions(suppress=True)

# Read date format from csv
def readDate(dfile):
    input_dates = list()
    julian_dates = list()
    year = 0
    with open('dateFormat.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for in_ in csv_reader:
            if line_count == 0:
                input_dates.append(in_)
            elif line_count == 1:
                julian_dates.append(in_)
            elif line_count == 2:
                year = in_
    return input_dates, julian_dates, year


## Setting file path parameters   _________________________________________________________________:
try:
    #print(sys.argv[1]) #tst
    sys.argv[1]
    if sys.argv.__contains__("-nogui"):
        print("\n\tAssuming file paths to be predefined...\n\t(\33[93mif this in not the case remove <-nogui> and run again\33[0m)")
        file_path_in = r"G:\My Drive\Stellenbosch\2022\716\ET\modeling\Data\in_"
        file_path_out = r"G:\My Drive\Stellenbosch\2022\716\ET\modeling\Data\out_"
        # file_path_in = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data"
        # file_path_out = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\output"
        rDate = readDate("dateFormat.csv")
        input_dates = rDate[0][0]
        julian_dates = rDate[0][1]
        year = rDate[0][2]
        print("\n\n", input_dates, "\n", julian_dates, "\n", year) #tst
    
    if sys.argv.__contains__("-lst"): # Not implemented
        print("\n\tUsing path list...")
    
    if sys.argv.__contains__("-mklst"): # Not implemented
        print("\n\tCreate new path list...")
except:
    #print("GUI") #tst
    root = tk.Tk()
    root.withdraw()

    file_path_in = filedialog.askdirectory(title="Please Select Input Folder")
    file_path_out = filedialog.askdirectory(title="Please Select Output Folder")
    input_dates = readDate("dateFormat.csv")
    print("\n\n input dates \n\n")
    print(input_dates)
    print("\n\n")
    
    
print(file_path_in, file_path_out, input_dates) #tst

## Clipping rasters to shapefile extent   _________________________________________________________:
def clipRast(outName, inRast, ext, reCreate=False):
    '''
    (rasterOut, rasterToClip, extent, save over)
    [path, path, path, bool]
    '''
    if not reCreate:
        if not os.path.exists(outName):
            ras = gdal.Open(inRast)
            srs = ras.GetProjectionRef()
            try:
                gdal.Warp(outName, ras, cutlineDSName=ext, cropToCutline=ras, srcSRS=srs)
                return [True, "Clipped Successfully"]
            except:
                return [False, "Could not Clip..."]
        else:
            return [True, "Clip Exists"]
    else:
        ras = gdal.Open(inRast)
        srs = ras.GetProjectionRef()
        try:
            gdal.Warp(outName, ras, cutlineDSName=ext, cropToCutline=ras, srcSRS=srs)
            return [True, "Clipped Successfully"]
        except:
            return [False, "Could not Clip..."]


def main(date, jdate):

    # define input files ______________________________________________________________________________:
    par = parm.PARAMS(file_path_in, file_path_out, date) #[0] <-should be able to handle any range a.t.m
    # clipPar = cparm.PARAMclip()
    print(par.getFilePathIN("time")) #tst

    # clip inputs files  ______________________________________________________________________________:
    print(par.getClipPathIN("albedo"), par.getFilePathIN("albedo"), par.getExtent()) #tst
    for i in par.getKeysIN():
        print(f'{i : >14}', end=' - ')
        state = clipRast(par.getClipPathIN(i), par.getFilePathIN(i), par.getExtent())
        if state[0]:
            print(f'{"done" : <7} : {state[1]}')
        else:
            print(f'{"failed" : <7} : {state[1]}')

    # read inputs files  ______________________________________________________________________________:
    dest_lst = gdal.Open(par.getClipPathIN("lst"))
    lst = dest_lst.GetRasterBand(1)
    nD = lst.GetNoDataValue()
    lst = lst.ReadAsArray()
    lst[lst == nD] = np.nan
    # print("LST ", lst, "\nmin: ", np.nanmin(lst), "\nmax: ", np.nanmax(lst)) #tst

    dest_albedo = gdal.Open(par.getClipPathIN("albedo"))
    r0 = dest_albedo.GetRasterBand(1)
    nD = r0.GetNoDataValue()
    r0 = r0.ReadAsArray()
    r0[r0 == nD] = np.nan
    # print("R0 ", r0, "\nmin: ", np.nanmin(r0), "\nmax: ", np.nanmax(r0)) #tst

    dest_ndvi = gdal.Open(par.getClipPathIN("ndvi"))
    ndvi = dest_ndvi.GetRasterBand(1)
    nD = ndvi.GetNoDataValue()
    ndvi = ndvi.ReadAsArray()
    ndvi[ndvi == nD] = np.nan
    # print("NDVI ", ndvi, "\nmin: ", np.nanmin(ndvi), "\nmax: ", np.nanmax(ndvi)) #tst

    desttime = gdal.Open(par.getClipPathIN("time"))
    dtime = desttime.GetRasterBand(1)
    nD = dtime.GetNoDataValue()
    dtime = dtime.ReadAsArray()
    dtime[dtime == nD] = np.nan
    # print("dtime ", dtime, "\nmin: ", np.nanmin(dtime), "\nmax: ", np.nanmax(dtime)) #tst

    dest_lat = gdal.Open(par.getClipPathIN("lat"))
    lat_deg = dest_lat.GetRasterBand(1)
    nD = lat_deg.GetNoDataValue()
    lat_deg = lat_deg.ReadAsArray()
    lat_deg[lat_deg == nD] = np.nan
    # print("Lat ", lat_deg, "\nmin: ", np.nanmin(lat_deg), "\nmax: ", np.nanmax(lat_deg)) #tst

    dest_lon = gdal.Open(par.getClipPathIN("lon"))
    lon_deg = dest_lon.GetRasterBand(1)
    nD = lon_deg.GetNoDataValue()
    lon_deg = lon_deg.ReadAsArray()
    lon_deg[lon_deg == nD] = np.nan
    # print("Lon ", lon_deg, "\nmin: ", np.nanmin(lon_deg), "\nmax: ", np.nanmax(lon_deg)) #tst

    dest_dem = gdal.Open(par.getClipPathIN("dem"))
    z = dest_dem.GetRasterBand(1)
    nD = z.GetNoDataValue()
    z = z.ReadAsArray()
    #z[np.isnan(lst)] = np.nan
    # z[z == nD] = np.nan
    # print("Z ", z, "\nmin: ", np.nanmin(z), "\nmax: ", np.nanmax(z)) #tst

    dest_slope = gdal.Open(par.getClipPathIN("slope"))
    slope_deg = dest_slope.GetRasterBand(1)
    nD = slope_deg.GetNoDataValue()
    slope_deg = slope_deg.ReadAsArray()
    slope_deg[slope_deg == nD] = np.nan
    # print("Slope ", slope_deg, "\nmin: ", np.nanmin(slope_deg), "\nmax: ", np.nanmax(slope_deg)) #tst

    dest_aspect = gdal.Open(par.getClipPathIN("aspect"))
    aspect_deg = dest_aspect.GetRasterBand(1)
    nD = aspect_deg.GetNoDataValue()
    aspect_deg = aspect_deg.ReadAsArray()
    aspect_deg[aspect_deg == nD] = np.nan
    # print("Aspect ", aspect_deg, "\nmin: ", np.nanmin(aspect_deg), "\nmax: ", np.nanmax(aspect_deg)) #tst

    dest_lm = gdal.Open(par.getClipPathIN("landMask"))
    land_mask = dest_lm.GetRasterBand(1)
    nD = land_mask.GetNoDataValue()
    land_mask = land_mask.ReadAsArray()
    # land_mask[np.isnan(lst)] = np.nan
    # land_mask[land_mask == nD] = np.nan
    # print("LandMask ", land_mask, "\nmin: ", np.nanmin(land_mask), "\nmax: ", np.nanmax(land_mask)) #tst

    #dest_bulk = gdal.Open(par.getClipPathIN("bulk"))
    #bulk = dest_bulk.GetRasterBand(1).ReadAsArray()

    dest_maxobs = gdal.Open(par.getClipPathIN("maxObs"))
    z_obst_max = dest_maxobs.GetRasterBand(1)
    nD = z_obst_max.GetNoDataValue()
    z_obst_max = z_obst_max.ReadAsArray()
    #z_obst_max[np.isnan(lst)] = np.nan
    # z_obst_max[z_obst_max == nD] = np.nan
    # print("ZobsMax ", z_obst_max, "\nmin: ", np.nanmin(z_obst_max), "\nmax: ", np.nanmax(z_obst_max)) #tst

    """ dest_pairsea24 = gdal.Open(par.getClipPathIN("pair_24_0"))
    p_air_0_24 = dest_pairsea24.GetRasterBand(1).ReadAsArray()
    p_air_0_24 = meteo.air_pressure_kpa2mbar(p_air_0_24)
    p_air_0_24[np.isnan(lst)] = np.nan

    dest_pairseainst = gdal.Open(par.getClipPathIN("pair_inst_0"))
    p_air_0_i = dest_pairseainst.GetRasterBand(1).ReadAsArray()
    p_air_0_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_i)
    p_air_0_i[np.isnan(lst)] = np.nan

    dest_pairinst = gdal.Open(par.getClipPathIN("pair_inst"))
    p_air_i = dest_pairinst.GetRasterBand(1).ReadAsArray()
    p_air_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_i)
    p_air_i[np.isnan(lst)] = np.nan

    dest_precip = gdal.Open(par.getClipPathIN("pre"))
    P_24 = dest_precip.GetRasterBand(1).ReadAsArray()
    P_24[np.isnan(lst)] = np.nan

    dest_hum24 = gdal.Open(par.getClipPathIN("hum_24"))
    qv_24 = dest_hum24.GetRasterBand(1).ReadAsArray()
    # qv_24[np.isnan(lst)] = np.nan

    dest_huminst = gdal.Open(par.getClipPathIN("hum_inst"))
    qv_i = dest_huminst.GetRasterBand(1).ReadAsArray()
    qv_i[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(par.getClipPathIN("tair_24"))
    t_air_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    #t_air_24 = ETLook.meteo.disaggregate_air_temperature_daily(t_air_24_coarse, z, z_coarse, lapse)
    # t_air_24[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(par.getClipPathIN("tair_max_24"))
    t_air_max_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    t_air_max_24[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(par.getClipPathIN("tair_min_24"))
    t_air_min_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    t_air_min_24[np.isnan(lst)] = np.nan

    dest_tairinst = gdal.Open(par.getClipPathIN("tair_inst"))
    t_air_i = dest_tairinst.GetRasterBand(1).ReadAsArray()
    # t_air_i[np.isnan(lst)] = np.nan

    dest_tairamp = gdal.Open(par.getClipPathIN("tair_amp"))
    t_amp_year = dest_tairamp.GetRasterBand(1).ReadAsArray()
    t_amp_year[np.isnan(lst)] = np.nan
    """
    dest_wind24 = gdal.Open(par.getClipPathIN("wind_24"))
    u_24 = dest_wind24.GetRasterBand(1)
    nD = u_24.GetNoDataValue()
    u_24 = u_24.ReadAsArray()
    # u_24[np.isnan(lst)] = np.nan
    # u_24[u_24 == nD] = np.nan

    dest_windinst = gdal.Open(par.getClipPathIN("wind_inst"))
    u_i = dest_windinst.GetRasterBand(1)
    nD = u_i.GetNoDataValue()
    u_i = u_i.ReadAsArray()
    # u_i[np.isnan(lst)] = np.nan
    # u_24[u_24 == nD] = np.nan

    dest_watcol = gdal.Open(par.getClipPathIN("watCol_inst"))
    wv_i = dest_watcol.GetRasterBand(1)
    nD = wv_i.GetNoDataValue()
    wv_i = wv_i.ReadAsArray()
    # wv_i[np.isnan(lst)] = np.nan
    # u_24[u_24 == nD] = np.nan

    dest_trans = gdal.Open(par.getClipPathIN("trans_24"))
    trans_24 = dest_trans.GetRasterBand(1)
    nD = trans_24.GetNoDataValue()
    trans_24 = trans_24.ReadAsArray()
    #trans_24[np.isnan(lst)] = np.nan
    # u_24[u_24 == nD] = np.nan
    # print("Trans24 ", trans_24, "\nmin: ", np.nanmin(trans_24), "\nmax: ", np.nanmax(trans_24)) #tst

	# define output filepaths
    # par_out = parm.PARAMS(file_path_out, date)

    # define prediction extent  ______________________________________________________________________:
    geo_ex = dest_lst.GetGeoTransform()
    proj_ex = dest_lst.GetProjection()

    # Create QC array  _______________________________________________________________________________:
    QC = np.ones(lst.shape)
    # QC[np.isnan(lst)] = np.nan


    # Parameters _________________________________________________________________Not Yet Adjustable__:

    doy = int(jdate)
    aod550_i = 0.01 # https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MOD04_L2 heb niet echt een standaard product hiervan gevonden
    se_top = 0.5
    porosity = 0.4

    # **effective_leaf_area_index
    # constants or predefined:
    vc_pow = 0.7
    lai_pow = -0.45

    # **atmospheric canopy resistance
    # constants or predefined:
    diffusion_slope = -1.33
    diffusion_intercept = 1.15
    t_opt = 25 # optimal temperature for plant growth
    t_min = 0 # minimal temperature for plant growth
    t_max = 50 # maximal temperature for plant growth
    vpd_slope = -0.3
    rs_min = 70
    rcan_max = 1000000

    # **net radiation canopy
    # constants or predefined:
    vp_slope = 0.14
    vp_offset = 0.34
    lw_slope = 1.35
    lw_offset = 0.35
    int_max = 0.2

    # **canopy resistance
    # constants or predefined:
    z_obs = 2
    z_b = 100
    z0m_bare = 0.001
    r0_bare = 0.38
    r0_full = 0.18
    tenacity = 1.5
    disp_bare = 0.0
    disp_full = 0.667
    fraction_h_bare = 0.65
    fraction_h_full = 0.95
    z0m_full = 0.1

    # **initial canopy aerodynamic resistance
    # constants or predefined:
    ndvi_obs_min = 0.25
    ndvi_obs_max = 0.75
    obs_fr = 0.25
    dem_resolution = 250

    # **ETLook.unstable.initial_friction_velocity_daily
    # constants or predefined:
    c1 = 1

    # **ETLook.unstable.transpiration
    # constants or predefined:
    iter_h = 3

    # **ETLook.resistance.soil_resistance
    # constants or predefined:
    r_soil_pow = -2.1
    r_soil_min = 800


    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily
    # constants or predefined:
    #porosity = 0.4 #Note: soil dependent
    #se_top = 1.0 #Note should be input !
    rn_slope = 0.92
    rn_offset = -61.0

    # **ETLook.unstable.evaporation
    # constants or predefined:
    r0_grass = 0.23

    #=================================================================================================:
    # LAI

    nd_min = np.nanmin(ndvi)
    nd_max = np.nanmax(ndvi)

    vc = leaf.vegetation_cover(ndvi, nd_min, nd_max, vc_pow)
    vc_min = np.nanmin(vc)
    vc_max = np.nanmax(vc)

    lai = leaf.leaf_area_index(vc, vc_min, vc_max, lai_pow)
    lai_eff = leaf.effective_leaf_area_index(lai)

#==========================================================================================================================

    # vc[np.isnan(QC)] = np.nan
    # lai[np.isnan(QC)] = np.nan
    # lai_eff[np.isnan(QC)] = np.nan
    if out.vc == 1:
        PF.Save_as_tiff(par.getFilePathOUT("vc"), vc, geo_ex, proj_ex)
    if out.lai == 1:
        PF.Save_as_tiff(par.getFilePathOUT("lai"), lai, geo_ex, proj_ex)
    if out.lai_eff == 1:
        PF.Save_as_tiff(par.getFilePathOUT("lai_eff"), lai_eff, geo_ex, proj_ex)


    #*******TRANSPIRATION COMPONENT****************************************************************

    # **soil fraction**************************************************************
    sf_soil = radiation.soil_fraction(lai)

    # sf_soil[np.isnan(QC)] = np.nan
    if out.sf_soil == 1:
        PF.Save_as_tiff(par.getFilePathOUT("sf_soil"), sf_soil, geo_ex, proj_ex)

    # **atmospheric canopy resistance***********************************************
    iesd = solar_radiation.inverse_earth_sun_distance(doy)
    sc = solar_radiation.seasonal_correction(doy)
    day_angle = clear_sky_radiation.day_angle(doy)
    decl = solar_radiation.declination(doy)
    lat = solar_radiation.latitude_rad(lat_deg)
    slope = solar_radiation.slope_rad(slope_deg)
    aspect = solar_radiation.aspect_rad(aspect_deg)

    print(sc)
    print(decl)
    print(iesd)
    print(lat)
    print(slope)
    print(aspect)

    print("\n\ncalc ra_24_toa\n\n")
    # print("sc: " + str(type(sc)))
    print("decl: " + str(type(decl)))
    print("\tndim: ", str(decl.ndim))
    print("\tshape: ", str(decl.shape))
    # print("iesd: " + str(type(iesd)))
    print("lat: " + str(type(lat)))
    print("\tndim: ", str(lat.ndim))
    print("\tshape: ", str(lat.shape))
    
    print("slope: " + str(type(slope)))
    print("\tndim: ", str(slope.ndim))
    print("\tshape: ", str(slope.shape))
    # print("aspect: " + str(type(aspect)))

    ra_24_toa = solar_radiation.daily_solar_radiation_toa(sc, decl, iesd, lat, slope, aspect)
    ws = solar_radiation.sunset_hour_angle(lat, decl)
    ra_24_toa_flat = solar_radiation.daily_solar_radiation_toa_flat(decl, iesd, lat, ws)
    diffusion_index = solar_radiation.diffusion_index(trans_24, diffusion_slope, diffusion_intercept)

    # choose one of the two options below
    #ra_24 = ETLook.solar_radiation.daily_solar_radiation_flat(ra_24_toa_flat, trans_24)
    ra_24 = solar_radiation.daily_total_solar_radiation(ra_24_toa, ra_24_toa_flat, diffusion_index, trans_24)
    stress_rad = stress.stress_radiation(ra_24)
    p_air_24 = meteo.air_pressure_daily(z, p_air_0_24)
    vp_24 = meteo.vapour_pressure_from_specific_humidity_daily(qv_24, p_air_24)
    svp_24 = meteo.saturated_vapour_pressure_average(
                meteo.saturated_vapour_pressure_maximum(t_air_max_24),
                meteo.saturated_vapour_pressure_minimum(t_air_min_24))
    vpd_24 = meteo.vapour_pressure_deficit_daily(svp_24, vp_24)
    stress_vpd = stress.stress_vpd(vpd_24, vpd_slope)
    stress_temp = stress.stress_temperature(t_air_24, t_opt, t_min, t_max)
    r_canopy_0 = resistance.atmospheric_canopy_resistance(lai_eff, stress_rad, stress_vpd, stress_temp, rs_min, rcan_max)

    # Save as tiff files
    # lat[np.isnan(QC)] = np.nan
    # slope[np.isnan(QC)] = np.nan
    # aspect[np.isnan(QC)] = np.nan
    # ra_24_toa[np.isnan(QC)] = np.nan
    # ws[np.isnan(QC)] = np.nan
    # ra_24_toa_flat[np.isnan(QC)] = np.nan
    # diffusion_index[np.isnan(QC)] = np.nan
    # ra_24[np.isnan(QC)] = np.nan
    # stress_rad[np.isnan(QC)] = np.nan
    # p_air_24[np.isnan(QC)] = np.nan
    # vp_24[np.isnan(QC)] = np.nan
    # svp_24[np.isnan(QC)] = np.nan
    # vpd_24[np.isnan(QC)] = np.nan
    # stress_vpd[np.isnan(QC)] = np.nan
    # stress_temp[np.isnan(QC)] = np.nan
    # r_canopy_0[np.isnan(QC)] = np.nan
    if out.lat == 1:
        PF.Save_as_tiff(par.getFilePathOUT("lat"), lat, geo_ex, proj_ex)
    if out.slope == 1:
        PF.Save_as_tiff(par.getFilePathOUT("slope"), slope, geo_ex, proj_ex)
    if out.aspect == 1:
        PF.Save_as_tiff(par.getFilePathOUT("aspect"), aspect, geo_ex, proj_ex)
    if out.ws == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ws"), ws, geo_ex, proj_ex)
    if out.ra_24_toa == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ra_24_toa"), ra_24_toa, geo_ex, proj_ex)
    if out.diffusion_index == 1:
        PF.Save_as_tiff(par.getFilePathOUT("diffusion_index"), diffusion_index, geo_ex, proj_ex)
    if out.ra_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ra_24"), ra_24, geo_ex, proj_ex)
    if out.stress_rad == 1:
        PF.Save_as_tiff(par.getFilePathOUT("stress_rad"), stress_rad, geo_ex, proj_ex)
    if out.p_air_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("p_air_24"), p_air_24, geo_ex, proj_ex)
    if out.vp_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("vp_24"), vp_24, geo_ex, proj_ex)
    if out.svp_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("svp_24"), svp_24, geo_ex, proj_ex)
    if out.vpd_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("vpd_24"), vpd_24, geo_ex, proj_ex)
    if out.stress_vpd == 1:
        PF.Save_as_tiff(par.getFilePathOUT("stress_vpd"), stress_vpd, geo_ex, proj_ex)
    if out.stress_temp == 1:
        PF.Save_as_tiff(par.getFilePathOUT("stress_temp"), stress_temp, geo_ex, proj_ex)
    if out.r_canopy_0 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("r_canopy_0"), r_canopy_0, geo_ex, proj_ex)

    # **net radiation canopy******************************************************
    t_air_k_24 = meteo.air_temperature_kelvin_daily(t_air_24)
    # select one of the below two
    #l_net = ETLook.radiation.longwave_radiation_fao_etref(t_air_k_24, vp_24, trans_24)
    l_net = radiation.longwave_radiation_fao(t_air_k_24, vp_24, trans_24, vp_slope, vp_offset, lw_slope, lw_offset)
    int_mm = evapotranspiration.interception_mm(P_24, vc, lai, int_max)
    lh_24 = meteo.latent_heat_daily(t_air_24)
    int_wm2 = radiation.interception_wm2(int_mm, lh_24)
    rn_24 = radiation.net_radiation(r0, ra_24, l_net, int_wm2)
    rn_24_canopy = radiation.net_radiation_canopy(rn_24, sf_soil)

    # Save as tiff files
    # t_air_k_24[np.isnan(QC)] = np.nan
    # l_net[np.isnan(QC)] = np.nan
    # int_mm[np.isnan(QC)] = np.nan
    # lh_24[np.isnan(QC)] = np.nan
    # int_wm2[np.isnan(QC)] = np.nan
    # rn_24[np.isnan(QC)] = np.nan
    # rn_24_canopy[np.isnan(QC)] = np.nan
    if out.t_air_k_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_air_k_24"), t_air_k_24, geo_ex, proj_ex)
    if out.l_net == 1:
        PF.Save_as_tiff(par.getFilePathOUT("l_net"), l_net, geo_ex, proj_ex)
    if out.int_mm == 1:
        PF.Save_as_tiff(par.getFilePathOUT("int_mm"), int_mm, geo_ex, proj_ex)
    if out.lh_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("lh_24"), lh_24, geo_ex, proj_ex)
    if out.int_wm2 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("int_wm2"), int_wm2, geo_ex, proj_ex)
    if out.rn_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rn_24"), rn_24, geo_ex, proj_ex)
    if out.rn_24_canopy == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rn_24_canopy"), rn_24_canopy, geo_ex, proj_ex)

    # **canopy resistance***********************************************************

    t_air_k_i = meteo.air_temperature_kelvin_inst(t_air_i)
    vp_i = meteo.vapour_pressure_from_specific_humidity_inst(qv_i, p_air_i)
    ad_moist_i = meteo.moist_air_density_inst(vp_i, t_air_k_i)
    ad_dry_i = meteo.dry_air_density_inst(p_air_i, vp_i, t_air_k_i)
    ad_i = meteo.air_density_inst(ad_dry_i, ad_moist_i)
    u_b_i_bare = soil_moisture.wind_speed_blending_height_bare(u_i, z0m_bare, z_obs, z_b)
    lon = solar_radiation.longitude_rad(lon_deg)
    ha = solar_radiation.hour_angle(sc, dtime, lon)
    I0 = clear_sky_radiation.solar_constant()
    ied = clear_sky_radiation.inverse_earth_sun_distance(day_angle)
    h0 = clear_sky_radiation.solar_elevation_angle(lat, decl, ha)
    h0ref = clear_sky_radiation.solar_elevation_angle_refracted(h0)
    m = clear_sky_radiation.relative_optical_airmass(p_air_i, p_air_0_i, h0ref)
    rotm = clear_sky_radiation.rayleigh_optical_thickness(m)
    Tl2 = clear_sky_radiation.linke_turbidity(wv_i, aod550_i, p_air_i, p_air_0_i)
    G0 = clear_sky_radiation.extraterrestrial_irradiance_normal(I0, ied)
    B0c = clear_sky_radiation.beam_irradiance_normal_clear(G0, Tl2, m, rotm, h0)
    Bhc = clear_sky_radiation.beam_irradiance_horizontal_clear(B0c, h0)
    Dhc = clear_sky_radiation.diffuse_irradiance_horizontal_clear(G0, Tl2, h0)
    ra_hor_clear_i = clear_sky_radiation.ra_clear_horizontal(Bhc, Dhc)
    emiss_atm_i = soil_moisture.atmospheric_emissivity_inst(vp_i, t_air_k_i)
    rn_bare = soil_moisture.net_radiation_bare(ra_hor_clear_i, emiss_atm_i, t_air_k_i, lst, r0_bare)
    rn_full = soil_moisture.net_radiation_full(ra_hor_clear_i, emiss_atm_i, t_air_k_i, lst, r0_full)
    h_bare = soil_moisture.sensible_heat_flux_bare(rn_bare, fraction_h_bare)
    h_full = soil_moisture.sensible_heat_flux_full(rn_full, fraction_h_full)
    u_b_i_full = soil_moisture.wind_speed_blending_height_full_inst(u_i, z0m_full, z_obs, z_b)
    u_star_i_bare = soil_moisture.friction_velocity_bare_inst(u_b_i_bare, z0m_bare, disp_bare, z_b)
    u_star_i_full = soil_moisture.friction_velocity_full_inst(u_b_i_full, z0m_full, disp_full, z_b)
    L_bare = soil_moisture.monin_obukhov_length_bare(h_bare, ad_i, u_star_i_bare, t_air_k_i)
    L_full = soil_moisture.monin_obukhov_length_full(h_full, ad_i, u_star_i_full, t_air_k_i)
    u_i_soil = soil_moisture.wind_speed_soil_inst(u_i, L_bare, z_obs)
    ras = soil_moisture.aerodynamical_resistance_soil(u_i_soil)
    raa = soil_moisture.aerodynamical_resistance_bare(u_i, L_bare, z0m_bare, disp_bare, z_obs)
    rac = soil_moisture.aerodynamical_resistance_full(u_i, L_full, z0m_full, disp_full, z_obs)
    t_max_bare = soil_moisture.maximum_temperature_bare(ra_hor_clear_i, emiss_atm_i, t_air_k_i, ad_i, raa, ras, r0_bare)
    t_max_full = soil_moisture.maximum_temperature_full(ra_hor_clear_i, emiss_atm_i, t_air_k_i, ad_i, rac, r0_full)
    w_i = soil_moisture.dew_point_temperature_inst(vp_i)
    t_dew_i = soil_moisture.dew_point_temperature_inst(vp_i)
    t_wet_i = soil_moisture.wet_bulb_temperature_inst(t_air_i, t_dew_i)
    t_wet_k_i = meteo.wet_bulb_temperature_kelvin_inst(t_wet_i)
    lst_max = soil_moisture.maximum_temperature(t_max_bare, t_max_full, vc)
    lst_min = soil_moisture.minimum_temperature(t_wet_k_i, t_air_k_i, vc)
    se_root = soil_moisture.soil_moisture_from_maximum_temperature(lst_max, lst, lst_min)
    stress_moist = stress.stress_moisture(se_root, tenacity)
    r_canopy_0 = resistance.atmospheric_canopy_resistance(lai_eff, stress_rad, stress_vpd, stress_temp, rs_min, rcan_max)
    r_canopy = resistance.canopy_resistance(r_canopy_0, stress_moist, rcan_max)

    # Save as tiff files
    # t_air_k_i[np.isnan(QC)] = np.nan
    # vp_i[np.isnan(QC)] = np.nan
    # ad_moist_i[np.isnan(QC)] = np.nan
    # ad_dry_i[np.isnan(QC)] = np.nan
    # ad_i[np.isnan(QC)] = np.nan
    # u_b_i_bare[np.isnan(QC)] = np.nan
    # lon[np.isnan(QC)] = np.nan
    # ha[np.isnan(QC)] = np.nan
    # h0[np.isnan(QC)] = np.nan
    # h0ref[np.isnan(QC)] = np.nan
    # m[np.isnan(QC)] = np.nan
    # rotm[np.isnan(QC)] = np.nan
    # Tl2[np.isnan(QC)] = np.nan
    # B0c[np.isnan(QC)] = np.nan
    # Bhc[np.isnan(QC)] = np.nan
    # Dhc[np.isnan(QC)] = np.nan
    # ra_hor_clear_i[np.isnan(QC)] = np.nan
    # emiss_atm_i[np.isnan(QC)] = np.nan
    # rn_bare[np.isnan(QC)] = np.nan
    # rn_full[np.isnan(QC)] = np.nan
    # u_b_i_full[np.isnan(QC)] = np.nan
    # u_star_i_bare[np.isnan(QC)] = np.nan
    # u_star_i_full[np.isnan(QC)] = np.nan
    # u_i_soil[np.isnan(QC)] = np.nan
    # ras[np.isnan(QC)] = np.nan
    # raa[np.isnan(QC)] = np.nan
    # rac[np.isnan(QC)] = np.nan
    # t_max_bare[np.isnan(QC)] = np.nan
    # t_max_full[np.isnan(QC)] = np.nan
    # w_i[np.isnan(QC)] = np.nan
    # t_dew_i[np.isnan(QC)] = np.nan
    # t_wet_i[np.isnan(QC)] = np.nan
    # t_wet_k_i[np.isnan(QC)] = np.nan
    # lst_max[np.isnan(QC)] = np.nan
    # se_root[np.isnan(QC)] = np.nan
    # stress_moist[np.isnan(QC)] = np.nan
    # r_canopy_0[np.isnan(QC)] = np.nan
    # r_canopy[np.isnan(QC)] = np.nan
    if out.t_air_k_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_air_k_i"), t_air_k_i, geo_ex, proj_ex)
    if out.vp_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("vp_i"), vp_i, geo_ex, proj_ex)
    if out.ad_moist_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ad_moist_i"), ad_moist_i, geo_ex, proj_ex)
    if out.ad_dry_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ad_dry_i"), ad_dry_i, geo_ex, proj_ex)
    if out.ad_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ad_i"), ad_i, geo_ex, proj_ex)
    if out.u_b_i_bare == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_b_i_bare"), u_b_i_bare, geo_ex, proj_ex)
    if out.lon == 1:
        PF.Save_as_tiff(par.getFilePathOUT("lon"), lon, geo_ex, proj_ex)
    if out.ha == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ha"), ha, geo_ex, proj_ex)
    if out.ied == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ied"), ied, geo_ex, proj_ex)
    if out.h0 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("h0"), h0, geo_ex, proj_ex)
    if out.h0ref == 1:
        PF.Save_as_tiff(par.getFilePathOUT("h0ref"), h0ref, geo_ex, proj_ex)
    if out.m == 1:
        PF.Save_as_tiff(par.getFilePathOUT("m"), m, geo_ex, proj_ex)
    if out.rotm == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rotm"), rotm, geo_ex, proj_ex)
    if out.Tl2 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("Tl2"), Tl2, geo_ex, proj_ex)
    if out.B0c == 1:
        PF.Save_as_tiff(par.getFilePathOUT("B0c"), B0c, geo_ex, proj_ex)
    if out.Bhc == 1:
        PF.Save_as_tiff(par.getFilePathOUT("Bhc"), Bhc, geo_ex, proj_ex)
    if out.Dhc == 1:
        PF.Save_as_tiff(par.getFilePathOUT("Dhc"), Dhc, geo_ex, proj_ex)
    if out.ra_hor_clear_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ra_hor_clear_i"), ra_hor_clear_i, geo_ex, proj_ex)
    if out.emiss_atm_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("emiss_atm_i"), emiss_atm_i, geo_ex, proj_ex)
    if out.rn_bare == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rn_bare"), rn_bare, geo_ex, proj_ex)
    if out.rn_full == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rn_full"), rn_full, geo_ex, proj_ex)
    if out.u_b_i_full == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_b_i_full"), u_b_i_full, geo_ex, proj_ex)
    if out.u_star_i_bare == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_star_i_bare"), u_star_i_bare, geo_ex, proj_ex)
    if out.u_star_i_full == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_star_i_full"), u_star_i_full, geo_ex, proj_ex)
    if out.u_i_soil == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_i_soil"), u_i_soil, geo_ex, proj_ex)
    if out.ras == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ras"), ras, geo_ex, proj_ex)
    if out.raa == 1:
        PF.Save_as_tiff(par.getFilePathOUT("raa"), raa, geo_ex, proj_ex)
    if out.rac == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rac"), rac, geo_ex, proj_ex)
    if out.t_max_bare == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_max_bare"), t_max_bare, geo_ex, proj_ex)
    if out.t_max_full == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_max_full"), t_max_full, geo_ex, proj_ex)
    if out.w_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("w_i"), w_i, geo_ex, proj_ex)
    if out.t_dew_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_dew_i"), t_dew_i, geo_ex, proj_ex)
    if out.t_wet_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_wet_i"), t_wet_i, geo_ex, proj_ex)
    if out.t_wet_k_i == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_wet_k_i"), t_wet_k_i, geo_ex, proj_ex)
    if out.lst_max == 1:
        PF.Save_as_tiff(par.getFilePathOUT("lst_max"), lst_max, geo_ex, proj_ex)
    if out.se_root == 1:
        PF.Save_as_tiff(par.getFilePathOUT("se_root"), se_root, geo_ex, proj_ex)
    if out.stress_moist == 1:
        PF.Save_as_tiff(par.getFilePathOUT("stress_moist"), stress_moist, geo_ex, proj_ex)
    if out.r_canopy_0 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("r_canopy_0"), r_canopy_0, geo_ex, proj_ex)
    if out.r_canopy == 1:
        PF.Save_as_tiff(par.getFilePathOUT("r_canopy"), r_canopy, geo_ex, proj_ex)

    # **initial canopy aerodynamic resistance***********************************************************

    z_obst = roughness.obstacle_height(ndvi, z_obst_max, ndvi_obs_min, ndvi_obs_max, obs_fr)
    z_oro = roughness.orographic_roughness(slope, dem_resolution) #careful - standard res is set to 250 # !!!
    z0m = roughness.roughness_length(lai, z_oro, z_obst, z_obst_max, land_mask)
    ra_canopy_init = neutral.initial_canopy_aerodynamic_resistance(u_24, z0m, z_obs)

    # Save as tiff files
    # z_obst[np.isnan(QC)] = np.nan
    # z_oro[np.isnan(QC)] = np.nan
    # z0m[np.isnan(QC)] = np.nan
    # ra_canopy_init[np.isnan(QC)] = np.nan
    if out.z_obst == 1:
        PF.Save_as_tiff(par.getFilePathOUT("z_obst"), z_obst, geo_ex, proj_ex)
    if out.z_oro == 1:
        PF.Save_as_tiff(par.getFilePathOUT("z_oro"), z_oro, geo_ex, proj_ex)
    if out.z0m == 1:
        PF.Save_as_tiff(par.getFilePathOUT("z0m"), z0m, geo_ex, proj_ex)
    if out.ra_canopy_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ra_canopy_init"), ra_canopy_init, geo_ex, proj_ex)

    # **windspeed blending height daily***********************************************************
    u_b_24 = meteo.wind_speed_blending_height_daily(u_24, z_obs, z_b)

    # Save as tiff files
    # u_b_24[np.isnan(QC)] = np.nan
    if out.u_b_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_b_24"), u_b_24, geo_ex, proj_ex)

    # **ETLook.unstable.initial_friction_velocity_daily***********************************************************
    disp = roughness.displacement_height(lai, z_obst, land_mask, c1)
    u_star_24_init = unstable.initial_friction_velocity_daily(u_b_24, z0m, disp, z_b)

    # Save as tiff files
    # disp[np.isnan(QC)] = np.nan
    # u_star_24_init[np.isnan(QC)] = np.nan
    if out.disp == 1:
        PF.Save_as_tiff(par.getFilePathOUT("disp"), disp, geo_ex, proj_ex)
    if out.u_star_24_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_star_24_init"), u_star_24_init, geo_ex, proj_ex)

    # **ETLook.neutral.initial_daily_transpiration***********************************************************
    ad_dry_24 = meteo.dry_air_density_daily(p_air_24, vp_24, t_air_k_24)
    ad_moist_24 = meteo.moist_air_density_daily(vp_24, t_air_k_24)
    ad_24 = meteo.air_density_daily(ad_dry_24, ad_moist_24)
    psy_24 = meteo.psychrometric_constant_daily(p_air_24, lh_24)
    ssvp_24 = meteo.slope_saturated_vapour_pressure_daily(t_air_24)
    t_24_init = neutral.initial_daily_transpiration(rn_24_canopy, ssvp_24, ad_24, vpd_24, psy_24, r_canopy, ra_canopy_init)

    # Save as tiff files
    # ad_dry_24[np.isnan(QC)] = np.nan
    # ad_moist_24[np.isnan(QC)] = np.nan
    # ad_24[np.isnan(QC)] = np.nan
    # psy_24[np.isnan(QC)] = np.nan
    # ssvp_24[np.isnan(QC)] = np.nan
    # t_24_init[np.isnan(QC)] = np.nan
    if out.ad_dry_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ad_dry_24"), ad_dry_24, geo_ex, proj_ex)
    if out.ad_moist_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ad_moist_24"), ad_moist_24, geo_ex, proj_ex)
    if out.ad_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ad_24"), ad_24, geo_ex, proj_ex)
    if out.psy_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("psy_24"), psy_24, geo_ex, proj_ex)
    if out.ssvp_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ssvp_24"), ssvp_24, geo_ex, proj_ex)
    if out.t_24_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_24_init"), t_24_init, geo_ex, proj_ex)

    # **ETLook.unstable.initial_sensible_heat_flux_canopy_daily***********************************************************
    h_canopy_24_init = unstable.initial_sensible_heat_flux_canopy_daily(rn_24_canopy, t_24_init)

    # Save as tiff files
    # h_canopy_24_init[np.isnan(QC)] = np.nan
    if out.h_canopy_24_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("h_canopy_24_init"), h_canopy_24_init, geo_ex, proj_ex)

    # **ETLook.unstable.transpiration***********************************************************

    t_24 = unstable.transpiration(rn_24_canopy, ssvp_24, ad_24, vpd_24, psy_24, r_canopy, h_canopy_24_init, t_air_k_24, u_star_24_init, z0m, disp, u_b_24, z_obs, z_b, iter_h)
    t_24_mm = unstable.transpiration_mm(t_24, lh_24)

    # Save as tiff files
    # t_24[np.isnan(QC)] = np.nan
    # t_24_mm[np.isnan(QC)] = np.nan
    if out.t_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_24"), t_24, geo_ex, proj_ex)
    if out.t_24_mm == 1:
        PF.Save_as_tiff(par.getFilePathOUT("t_24_mm"), t_24_mm, geo_ex, proj_ex)

    #*******EVAPORATION COMPONENT****************************************************************

    # **ETLook.radiation.net_radiation_soil***********************************************************
    sf_soil = radiation.soil_fraction(lai)
    rn_24_soil = radiation.net_radiation_soil(rn_24, sf_soil)

    # Save as tiff files
    # sf_soil[np.isnan(QC)] = np.nan
    # rn_24_soil[np.isnan(QC)] = np.nan
    if out.sf_soil == 1:
        PF.Save_as_tiff(par.getFilePathOUT("sf_soil"), sf_soil, geo_ex, proj_ex)
    if out.rn_24_soil == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rn_24_soil"), rn_24_soil, geo_ex, proj_ex)

    # **ETLook.resistance.soil_resistance***********************************************************

    r_soil = resistance.soil_resistance(se_top, land_mask, r_soil_pow, r_soil_min)

    # Save as tiff files
    # r_soil[np.isnan(QC)] = np.nan
    if out.r_soil == 1:
        PF.Save_as_tiff(par.getFilePathOUT("r_soil"), r_soil, geo_ex, proj_ex)

    # **ETLook.resistance.soil_resistance***********************************************************

    ra_soil_init = neutral.initial_soil_aerodynamic_resistance(u_24, z_obs)

    # Save as tiff files
    # ra_soil_init[np.isnan(QC)] = np.nan
    if out.ra_soil_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("ra_soil_init"), ra_soil_init, geo_ex, proj_ex)

    # **ETLook.meteo.wind_speed_blending_height_daily***********************************************************

    u_b_24 = meteo.wind_speed_blending_height_daily(u_24, z_obs, z_b)

    # Save as tiff files
    # u_b_24[np.isnan(QC)] = np.nan
    if out.u_b_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_b_24"), u_b_24, geo_ex, proj_ex)

    # **ETLook.unstable.initial_friction_velocity_soil_daily***********************************************************

    u_star_24_soil_init = unstable.initial_friction_velocity_soil_daily(u_b_24, disp, z_b)

    # Save as tiff files
    # u_star_24_soil_init[np.isnan(QC)] = np.nan
    if out.u_star_24_soil_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("u_star_24_soil_init"), u_star_24_soil_init, geo_ex, proj_ex)

    # **ETLook.unstable.initial_sensible_heat_flux_soil_daily***********************************************************

    stc = radiation.soil_thermal_conductivity(se_top)
    vhc = radiation.volumetric_heat_capacity(se_top, porosity)
    dd = radiation.damping_depth(stc, vhc)
    g0_bs = radiation.bare_soil_heat_flux(doy, dd, stc, t_amp_year, lat)
    g0_24 = radiation.soil_heat_flux(g0_bs, sf_soil, land_mask, rn_24_soil, trans_24, ra_24, l_net, rn_slope, rn_offset)
    e_24_init = neutral.initial_daily_evaporation(rn_24_soil, g0_24, ssvp_24, ad_24, vpd_24, psy_24, r_soil, ra_soil_init)
    h_soil_24_init = unstable.initial_sensible_heat_flux_soil_daily(rn_24_soil, e_24_init, g0_24)

    # Save as tiff files
    # g0_bs[np.isnan(QC)] = np.nan
    # g0_24[np.isnan(QC)] = np.nan
    # e_24_init[np.isnan(QC)] = np.nan
    # h_soil_24_init[np.isnan(QC)] = np.nan
    if out.g0_bs == 1:
        PF.Save_as_tiff(par.getFilePathOUT("g0_bs"), g0_bs, geo_ex, proj_ex)
    if out.g0_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("g0_24"), g0_24, geo_ex, proj_ex)
    if out.e_24_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("e_24_init"), e_24_init, geo_ex, proj_ex)
    if out.h_soil_24_init == 1:
        PF.Save_as_tiff(par.getFilePathOUT("h_soil_24_init"), h_soil_24_init, geo_ex, proj_ex)

    # **ETLook.unstable.evaporation***********************************************************

    e_24 = unstable.evaporation(rn_24_soil, g0_24, ssvp_24, ad_24, vpd_24, psy_24, r_soil, h_soil_24_init, t_air_k_24, u_star_24_soil_init, disp, u_b_24, z_b, z_obs, iter_h)
    e_24_mm = unstable.evaporation_mm(e_24, lh_24)
    et_24_mm = evapotranspiration.et_actual_mm(e_24_mm, t_24_mm)

    # Save as tiff files
    # e_24[np.isnan(QC)] = np.nan
    # e_24_mm[np.isnan(QC)] = np.nan
    # et_24_mm[np.isnan(QC)] = np.nan
    if out.e_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("e_24"), e_24, geo_ex, proj_ex)
    if out.e_24_mm == 1:
        PF.Save_as_tiff(par.getFilePathOUT("e_24_mm"), e_24_mm, geo_ex, proj_ex)
    if out.et_24_mm == 1:
        PF.Save_as_tiff(par.getFilePathOUT("et_24_mm"), et_24_mm, geo_ex, proj_ex)

    # **ETLook.unstable.evaporation***********************************************************
    rn_24_grass = radiation.net_radiation_grass(ra_24, l_net, r0_grass)
    et_ref_24 = evapotranspiration.et_reference(rn_24_grass, ad_24, psy_24, vpd_24, ssvp_24, u_24)
    et_ref_24_mm = evapotranspiration.et_reference_mm(et_ref_24, lh_24)

    # Save as tiff files
    # rn_24_grass[np.isnan(QC)] = np.nan
    # et_ref_24[np.isnan(QC)] = np.nan
    # et_ref_24_mm[np.isnan(QC)] = np.nan
    if out.rn_24_grass == 1:
        PF.Save_as_tiff(par.getFilePathOUT("rn_24_grass"), rn_24_grass, geo_ex, proj_ex)
    if out.et_ref_24 == 1:
        PF.Save_as_tiff(par.getFilePathOUT("et_ref_24"), et_ref_24, geo_ex, proj_ex)
    if out.et_ref_24_mm == 1:
        PF.Save_as_tiff(par.getFilePathOUT("et_ref_24_mm"), et_ref_24_mm, geo_ex, proj_ex)



# for i in input_dates:
for i in range(0,1): # temp
    print("Currently processing {", input_dates[i],"} ...")
    main(input_dates[i], julian_dates[i])