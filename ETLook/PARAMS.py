from typing import Dict


class PARAMS:
    def __init__(self, folder_path, input_dates:None):

        import os

        self.parDict = {
            "albedo" : os.path.join(os.path.join(folder_path, "Albedo_outputs/"), "Albedo_%s.tif" %input_dates),
            "ndvi" : os.path.join(os.path.join(folder_path, "Indecies/"), "Indeces_%s.tif" %input_dates),
            "lst" : os.path.join(os.path.join(folder_path, "LST/"), "LST_%s.tif" %input_dates),
            #"time" : os.path.join(folder_path, "Time_%s.tif"),
            "lat" : os.path.join(os.path.join(folder_path, "Lat_Long_Rasters/"), "Lat.tif"),
            "lon" : os.path.join(os.path.join(folder_path, "Lat_Long_Rasters/"), "Lon.tif"),
            "dem" : os.path.join(os.path.join(folder_path, "DEM_Derviatives/"), "dem.tif"),
            "slope" : os.path.join(os.path.join(folder_path, "DEM_Derviatives/"), "slope.tif"),
            "aspect" : os.path.join(os.path.join(folder_path, "DEM_Derviatives/"), "aspect.tif"),
            "landMask" : os.path.join(os.path.join(folder_path, "Land_Mask/"), "LandMask.tif"),
            "bulk" : os.path.join(os.path.join(folder_path, "LULC_Products/"), "BSR_AOI.tif"),
            "maxObs" : os.path.join(os.path.join(folder_path, "LULC_Products/"), "Max_height_AOI.tif"),
            #"pair_24_0" : "Pair_24_0_%s.tif",
            #"pair_inst_0" : "Pair_inst_0_%s.tif",
            #"pair_inst" : "Pair_inst_%s.tif",
            #"pre" : "Precipitation_%s.tif",
            #"hum_24" : "qv_24_%s.tif",
            #"hum_inst" : "qv_inst_%s.tif",
            #"tair_24" : "tair_24_%s.tif",
            #"tair_max_24" : "tair_max_24_%s.tif",
            #"tair_min_24" : "tair_min_24_%s.tif",
            #"tair_inst" : "tair_inst_%s.tif",
            #"tair_amp" : "Tair_amp_%s.tif",
            #"wind_24" : "wind_24_%s.tif",
            #"wind_inst" : "wind_inst_%s.tif",
            #"watCol_inst" : "wv_inst_%s.tif",
            "trans_24" : os.path.join(os.path.join(folder_path, "transmissivity_false/"), "transmissivity.tif"),

			# outputs
			"vc" : os.path.join(folder_path, "vc_%s.tif"),
			"lai" : os.path.join(folder_path, "lai_%s.tif" %input_dates),
			"lai_eff" : os.path.join(folder_path, "lai_eff_%s.tif" %input_dates),
			"sf_soil" : os.path.join(folder_path, "sf_soil_%s.tif" %input_dates),
			"lat" : os.path.join(folder_path, "lat_%s.tif" %input_dates),
			"slope" : os.path.join(folder_path, "slope_%s.tif" %input_dates),
			"aspect" : os.path.join(folder_path, "aspect_%s.tif" %input_dates),
			"ra_24_toa" : os.path.join(folder_path, "ra_24_toa_%s.tif" %input_dates),
			"ws": os.path.join(folder_path, "ws_%s.tif" %input_dates),
			"diffusion_index" : os.path.join(folder_path, "diffusion_index_%s.tif" %input_dates),
			"ra_24" : os.path.join(folder_path, "ra_24_%s.tif" %input_dates),
			"stress_rad" : os.path.join(folder_path, "stress_rad_%s.tif" %input_dates),
			"p_air_24" : os.path.join(folder_path, "_%s.tif" %input_dates),
			"vp_24" : os.path.join(folder_path, "vp_24_%s.tif" %input_dates),
			"svp_24" : os.path.join(folder_path, "svp_24_%s.tif" %input_dates),
			"vpd_24" : os.path.join(folder_path, "vpd_24_%s.tif" %input_dates),
			"stress_vpd" : os.path.join(folder_path, "stress_vpd_%s.tif" %input_dates),
			"stress_temp" : os.path.join(folder_path, "stress_temp_%s.tif" %input_dates),
			"r_canopy_0" : os.path.join(folder_path, "r_canopy_0_%s.tif" %input_dates),
			"t_air_k_24" : os.path.join(folder_path, "t_air_k_24_%s.tif" %input_dates),
			"l_net" : os.path.join(folder_path, "l_net_%s.tif" %input_dates),
			"int_mm" : os.path.join(folder_path, "int_mm_%s.tif" %input_dates),
			"lh_24" : os.path.join(folder_path, "lh_24_%s.tif" %input_dates),
			"int_wm2" : os.path.join(folder_path, "int_wm2_%s.tif" %input_dates),
			"rn_24" : os.path.join(folder_path, "rn_24_%s.tif" %input_dates),
			"rn_24_canopy" : os.path.join(folder_path, "rn_24_canopy_%s.tif" %input_dates),
			"t_air_k_i" : os.path.join(folder_path, "t_air_k_i_%s.tif" %input_dates),
			"vp_i" : os.path.join(folder_path, "vp_i_%s.tif" %input_dates),
			"ad_moist_i" : os.path.join(folder_path, "ad_moist_i_%s.tif" %input_dates),
			"ad_dry_i" : os.path.join(folder_path, "ad_dry_i_%s.tif" %input_dates),
			"ad_i" : os.path.join(folder_path, "ad_i_%s.tif" %input_dates),
			"u_b_i_bare" : os.path.join(folder_path, "u_b_i_bare_%s.tif" %input_dates),
			"lon" : os.path.join(folder_path, "lon_%s.tif" %input_dates),
			"ha" : os.path.join(folder_path, "ha_%s.tif" %input_dates),
			"ied" : os.path.join(folder_path, "ied_%s.tif" %input_dates),
			"h0" : os.path.join(folder_path, "h0_%s.tif" %input_dates),
			"h0ref" : os.path.join(folder_path, "h0ref_%s.tif" %input_dates),
			"m" : os.path.join(folder_path, "m_%s.tif" %input_dates),
			"rotm" : os.path.join(folder_path, "rotm_%s.tif" %input_dates),
			"Tl2"  : os.path.join(folder_path, "Tl2_%s.tif" %input_dates),
			"B0c" : os.path.join(folder_path, "B0c_%s.tif" %input_dates),
			"Bhc" : os.path.join(folder_path, "Bhc_%s.tif" %input_dates),
			"Dhc" : os.path.join(folder_path, "Dhc_%s.tif" %input_dates),
			"ra_hor_clear_i" : os.path.join(folder_path, "ra_hor_clear_i_%s.tif" %input_dates),
			"emiss_atm_i"  : os.path.join(folder_path, "emiss_atm_i_%s.tif" %input_dates),
			"rn_bare" : os.path.join(folder_path, "rn_bare_%s.tif" %input_dates),
			"rn_full" : os.path.join(folder_path, "rn_full_%s.tif" %input_dates),
			"u_b_i_full" : os.path.join(folder_path, "u_b_i_full_%s.tif" %input_dates),
			"u_star_i_bare" : os.path.join(folder_path, "u_star_i_bare_%s.tif" %input_dates),
			"u_star_i_full" : os.path.join(folder_path, "u_star_i_full_%s.tif" %input_dates),
			"u_i_soil" : os.path.join(folder_path, "u_i_soil_%s.tif" %input_dates),
			"ras"  : os.path.join(folder_path, "ras_%s.tif" %input_dates),
			"raa" : os.path.join(folder_path, "raa_%s.tif" %input_dates),
			"rac" : os.path.join(folder_path, "rac_%s.tif" %input_dates),
			"t_max_bare" : os.path.join(folder_path, "t_max_bare_%s.tif" %input_dates),
			"t_max_full" : os.path.join(folder_path, "t_max_full_%s.tif" %input_dates),
			"w_i" : os.path.join(folder_path, "w_i_%s.tif" %input_dates),
			"t_dew_i" : os.path.join(folder_path, "t_dew_i_%s.tif" %input_dates),
			"t_wet_i" : os.path.join(folder_path, "t_wet_i_%s.tif" %input_dates),
			"t_wet_k_i" : os.path.join(folder_path, "t_wet_k_i_%s.tif" %input_dates),
			"lst_max"  : os.path.join(folder_path, "lst_max_%s.tif" %input_dates),
			"se_root" : os.path.join(folder_path, "se_root_%s.tif" %input_dates),
			"stress_moist" : os.path.join(folder_path, "stress_moist_%s.tif" %input_dates),
			"r_canopy_0" : os.path.join(folder_path, "r_canopy_0_%s.tif" %input_dates),
			"r_canopy" : os.path.join(folder_path, "r_canopy_%s.tif" %input_dates),
			"z_obst" : os.path.join(folder_path, "z_obst_%s.tif" %input_dates),
			"z_oro"  : os.path.join(folder_path, "z_oro_%s.tif" %input_dates),
			"z0m" : os.path.join(folder_path, "z0m_%s.tif" %input_dates),
			"ra_canopy_init" : os.path.join(folder_path, "ra_canopy_init_%s.tif" %input_dates),
			"u_b_24" : os.path.join(folder_path, "u_b_24_%s.tif" %input_dates),
			"disp" : os.path.join(folder_path, "disp_%s.tif" %input_dates),
			"u_star_24_init" : os.path.join(folder_path, "u_star_24_init_%s.tif" %input_dates),
			"ad_dry_24" : os.path.join(folder_path, "ad_dry_24_%s.tif" %input_dates),
			"ad_moist_24" : os.path.join(folder_path, "ad_moist_24_%s.tif" %input_dates),
			"ad_24" : os.path.join(folder_path, "ad_24_%s.tif" %input_dates),
			"psy_24" : os.path.join(folder_path, "psy_24_%s.tif" %input_dates),
			"ssvp_24" : os.path.join(folder_path, "ssvp_24_%s.tif" %input_dates),
			"t_24_init" : os.path.join(folder_path, "t_24_init_%s.tif" %input_dates),
			"h_canopy_24_init" : os.path.join(folder_path, "h_canopy_24_init_%s.tif" %input_dates),
			"t_24" : os.path.join(folder_path, "t_24_%s.tif" %input_dates),
			"t_24_mm" : os.path.join(folder_path, "t_24_mm_%s.tif" %input_dates),
			"sf_soil" : os.path.join(folder_path, "sf_soil_%s.tif" %input_dates),
			"rn_24_soil": os.path.join(folder_path, "rn_24_soil_%s.tif" %input_dates),
			"r_soil" : os.path.join(folder_path, "r_soil_%s.tif" %input_dates),
			"ra_soil_init" : os.path.join(folder_path, "ra_soil_init_%s.tif" %input_dates),
			"u_b_24" : os.path.join(folder_path, "u_b_24_%s.tif" %input_dates),
			"u_star_24_soil_init" : os.path.join(folder_path, "u_star_24_soil_init_%s.tif" %input_dates),
			"g0_bs" : os.path.join(folder_path, "g0_bs_%s.tif" %input_dates),
			"g0_24" : os.path.join(folder_path, "g0_24_%s.tif" %input_dates),
			"e_24_init" : os.path.join(folder_path, "e_24_init_%s.tif" %input_dates),
			"h_soil_24_init" : os.path.join(folder_path, "h_soil_24_init_%s.tif" %input_dates),
			"e_24" : os.path.join(folder_path, "e_24_%s.tif" %input_dates),
			"e_24_mm" : os.path.join(folder_path, "e_24_mm_%s.tif" %input_dates),
			"et_24_mm" : os.path.join(folder_path, "et_24_mm_%s.tif" %input_dates),
			"rn_24_grass" : os.path.join(folder_path, "rn_24_grass_%s.tif" %input_dates),
			"et_ref_24" : os.path.join(folder_path, "et_ref_24_%s.tif" %input_dates),
			"et_ref_24_mm" : os.path.join(folder_path, "et_ref_24_mm_%s.tif" %input_dates)
        }

    def getFileName(self, key):
        return self.parDict.get(key)
    
