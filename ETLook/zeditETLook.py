from cmath import inf
import os
import sys
# from osgeo import osr
from osgeo import gdal
import matplotlib.pyplot as plt
import numpy as np

import outputs as out
import PARAMS as parm
import tkinter as tk
from tkinter import filedialog
import csv
import math
import zeditETLook, solar_radiation, clear_sky_radiation, meteo, radiation, evapotranspiration, soil_moisture, leaf, stress, resistance, roughness, neutral, unstable, outputs

""" current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)
print(parent)
import ETLook as ET
import Functions as PF """

#np.set_printoptions(suppress=True)

# Read date format from csv
def readDate(dfile):
    input_dates = list()
    with open('dateFormat.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for input_dates in csv_reader:
            if line_count == 0:
                print(f'Date formats are {", ".join(input_dates)}')
    return input_dates

## Setting file path parameters   _________________________________________________________________:
try:
    #print(sys.argv[1]) #tst
    sys.argv[1]
    if sys.argv.__contains__("-nogui"):
        print("\n\tAssuming file paths to be predefined...\n\t(\33[93mif this in not the case remove <-nogui> and run again\33[0m)")
        file_path_in = r"G:\My Drive\Stellenbosch\2022\716\ET\modeling\Data\in_"
        file_path_out = r"G:\My Drive\Stellenbosch\2022\716\ET\modeling\Data\out_"
        input_dates = readDate("dateFormat.csv")
    
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
    
    
print(file_path_in, file_path_out, input_dates) #tst

def main(date):

    # define input files ______________________________________________________________________________:
    par = parm.PARAMS(file_path_in, date) #[0] <-- limited to one date a.t.m
    print(par.getFileName("albedo")) #tst


    # read inputs files  ______________________________________________________________________________:
    dest_lst = gdal.Open(par.getFileName("lst"))
    lst = dest_lst.GetRasterBand(1).ReadAsArray()
    lst = np.round(lst, 4)
    lst[lst == -inf] = np.nan
    print("LST ", lst, "\nmin: ", np.nanmin(lst), "\nmax: ", np.nanmax(lst)) #tst

    dest_albedo = gdal.Open(par.getFileName("albedo"))
    r0 = dest_albedo.GetRasterBand(1).ReadAsArray()
    r0 = np.round(r0, 4)
    r0[r0 == inf] = np.nan
    print("R0 ", r0, "\nmin: ", np.nanmin(r0), "\nmax: ", np.nanmax(r0)) #tst

    dest_ndvi = gdal.Open(par.getFileName("ndvi"))
    ndvi = dest_ndvi.GetRasterBand(1).ReadAsArray()
    ndvi = np.round(ndvi, 4)
    ndvi[ndvi == -inf] = np.nan
    print("NDVI ", ndvi, "\nmin: ", np.nanmin(ndvi), "\nmax: ", np.nanmax(ndvi)) #tst

    """ desttime = gdal.Open(par.getFileName("time"))
    dtime = desttime.GetRasterBand(1).ReadAsArray()
    dtime[np.isnan(lst)] = np.nan """

    dest_lat = gdal.Open(par.getFileName("lat"))
    lat_deg = dest_lat.GetRasterBand(1).ReadAsArray()
    print("Lat ", lat_deg, "\nmin: ", np.nanmin(lat_deg), "\nmax: ", np.nanmax(lat_deg))

    dest_lon = gdal.Open(par.getFileName("lon"))
    lon_deg = dest_lon.GetRasterBand(1).ReadAsArray()
    print("Lon ", lon_deg, "\nmin: ", np.nanmin(lon_deg), "\nmax: ", np.nanmax(lon_deg))

    dest_dem = gdal.Open(par.getFileName("dem"))
    z = dest_dem.GetRasterBand(1).ReadAsArray()
    #z[np.isnan(lst)] = np.nan
    print("Z ", z, "\nmin: ", np.nanmin(z), "\nmax: ", np.nanmax(z))

    dest_slope = gdal.Open(par.getFileName("slope"))
    slope_deg = dest_slope.GetRasterBand(1).ReadAsArray()
    slope_deg = np.round(slope_deg, 4)
    slope_deg[slope_deg == inf] = np.nan
    print("Slope ", slope_deg, "\nmin: ", np.nanmin(slope_deg), "\nmax: ", np.nanmax(slope_deg))

    dest_aspect = gdal.Open(par.getFileName("aspect"))
    aspect_deg = dest_aspect.GetRasterBand(1).ReadAsArray()
    aspect_deg = np.round(aspect_deg, 4)
    aspect_deg[aspect_deg == inf] = np.nan
    print("Aspect ", aspect_deg, "\nmin: ", np.nanmin(aspect_deg), "\nmax: ", np.nanmax(aspect_deg))

    dest_lm = gdal.Open(par.getFileName("landMask"))
    land_mask = dest_lm.GetRasterBand(1).ReadAsArray()
    #land_mask[np.isnan(lst)] = np.nan
    print("LandMask ", land_mask, "\nmin: ", np.nanmin(land_mask), "\nmax: ", np.nanmax(land_mask))

    #dest_bulk = gdal.Open(par.getFileName("bulk"))
    #bulk = dest_bulk.GetRasterBand(1).ReadAsArray()

    dest_maxobs = gdal.Open(par.getFileName("maxObs"))
    z_obst_max = dest_maxobs.GetRasterBand(1).ReadAsArray()
    #z_obst_max[np.isnan(lst)] = np.nan
    print("ZobsMax ", z_obst_max, "\nmin: ", np.nanmin(z_obst_max), "\nmax: ", np.nanmax(z_obst_max))

    """ dest_pairsea24 = gdal.Open(par.getFileName("pair_24_0"))
    p_air_0_24 = dest_pairsea24.GetRasterBand(1).ReadAsArray()
    p_air_0_24 = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_24)
    p_air_0_24[np.isnan(lst)] = np.nan

    dest_pairseainst = gdal.Open(par.getFileName("pair_inst_0"))
    p_air_0_i = dest_pairseainst.GetRasterBand(1).ReadAsArray()
    p_air_0_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_0_i)
    p_air_0_i[np.isnan(lst)] = np.nan

    dest_pairinst = gdal.Open(par.getFileName("pair_inst"))
    p_air_i = dest_pairinst.GetRasterBand(1).ReadAsArray()
    p_air_i = ETLook.meteo.air_pressure_kpa2mbar(p_air_i)
    p_air_i[np.isnan(lst)] = np.nan

    dest_precip = gdal.Open(par.getFileName("pre"))
    P_24 = dest_precip.GetRasterBand(1).ReadAsArray()
    P_24[np.isnan(lst)] = np.nan

    dest_hum24 = gdal.Open(par.getFileName("hum_24"))
    qv_24 = dest_hum24.GetRasterBand(1).ReadAsArray()
    qv_24[np.isnan(lst)] = np.nan

    dest_huminst = gdal.Open(par.getFileName("hum_inst"))
    qv_i = dest_huminst.GetRasterBand(1).ReadAsArray()
    qv_i[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(par.getFileName("tair_24"))
    t_air_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    #t_air_24 = ETLook.meteo.disaggregate_air_temperature_daily(t_air_24_coarse, z, z_coarse, lapse)
    t_air_24[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(par.getFileName("tair_max_24"))
    t_air_max_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    t_air_max_24[np.isnan(lst)] = np.nan

    dest_tair24 = gdal.Open(par.getFileName("tair_min_24"))
    t_air_min_24 = dest_tair24.GetRasterBand(1).ReadAsArray()
    t_air_min_24[np.isnan(lst)] = np.nan

    dest_tairinst = gdal.Open(par.getFileName("tair_inst"))
    t_air_i = dest_tairinst.GetRasterBand(1).ReadAsArray()
    t_air_i[np.isnan(lst)] = np.nan

    dest_tairamp = gdal.Open(par.getFileName("tair_amp"))
    t_amp_year = dest_tairamp.GetRasterBand(1).ReadAsArray()
    t_amp_year[np.isnan(lst)] = np.nan

    dest_wind24 = gdal.Open(par.getFileName("wind_24"))
    u_24 = dest_wind24.GetRasterBand(1).ReadAsArray()
    u_24[np.isnan(lst)] = np.nan

    dest_windinst = gdal.Open(par.getFileName("wind_inst"))
    u_i = dest_windinst.GetRasterBand(1).ReadAsArray()
    u_i[np.isnan(lst)] = np.nan

    dest_watcol = gdal.Open(par.getFileName("watCol_inst"))
    wv_i = dest_watcol.GetRasterBand(1).ReadAsArray()
    wv_i[np.isnan(lst)] = np.nan """

    dest_trans = gdal.Open(par.getFileName("trans_24"))
    trans_24 = dest_trans.GetRasterBand(1).ReadAsArray()
    #trans_24[np.isnan(lst)] = np.nan
    print("Trans24 ", trans_24, "\nmin: ", np.nanmin(trans_24), "\nmax: ", np.nanmax(trans_24))

	# define output filepaths
    par_out = parm.PARAMS(file_path_out, date)

    # define prediction extent  ______________________________________________________________________:
    geo_ex = dest_lst.GetGeoTransform()
    proj_ex = dest_lst.GetProjection()

    # Create QC array  _______________________________________________________________________________:
    QC = np.ones(lst.shape)
    QC[np.isnan(lst)] = np.nan


    # Parameters _________________________________________________________________Not Yet Adjustable__:

    # **effective_leaf_area_index
    # constants or predefined:
    nd_min = 0.125
    nd_max = 0.8
    vc_pow = 0.7
    vc_min = 0
    vc_max = 0.9677324224821418
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

    nd_min = np.nanmin(ndvi)
    nd_max = np.nanmax(ndvi)

    vc = leaf.vegetation_cover(ndvi, nd_min, nd_max, vc_pow)
    vc_min = np.nanmin(vc)
    vc_max = np.nanmax(vc)

    lai = leaf.leaf_area_index(vc, vc_min, vc_max, lai_pow)
    lai_eff = leaf.effective_leaf_area_index(lai)

    """ plt.imshow(lai_eff)
    plt.show() """

# for i in input_dates:
for i in range(0,1): # temp
    print("Currently processing {", input_dates[i],"} ...")
    main(input_dates[i])