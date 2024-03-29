class PARAMS:
	def __init__(self, input_path:str, output_path:str, input_dates:str):

		"""
		Setting up file paths for use in code, takes: 
		(path to input data directory: str, path to output data directory: str, optional: date range format: str)
		returns --> class
		"""

		import os

		self.extent =  os.path.join(os.path.join(input_path, "Extent/"), "extent.shp")

		self.parDictClip = {
		"albedo" : os.path.join(os.path.join(input_path, "Albedo/"), "albedo_clip_%s.tif" %input_dates),
		"ndvi" : os.path.join(os.path.join(input_path, "Indecies/"), "Indices_clip_%s.tif" %input_dates),
		"lst" : os.path.join(os.path.join(input_path, "LST/"), "lst_clip_%s.tif" %input_dates),
		"time" : os.path.join(os.path.join(input_path, "Time/"), "solarTime_clip.tif"),
		"lat" : os.path.join(os.path.join(input_path, "LatLon/"), "lat_clip.tif"),
		"lon" : os.path.join(os.path.join(input_path, "LatLon/"), "lon_clip.tif"),
		"dem" : os.path.join(os.path.join(input_path, "Z/"), "dem_clip.tif"),
		"slope" : os.path.join(os.path.join(input_path, "Z/"), "slope_clip.tif"),
		"aspect" : os.path.join(os.path.join(input_path, "Z/"), "aspect_clip.tif"),
		"landMask" : os.path.join(os.path.join(input_path, "LandMask/"), "landMask_clip.tif"),
		"bulk" : os.path.join(os.path.join(input_path, "LULC_Products/"), "BRS_clip.tif"),
		"maxObs" : os.path.join(os.path.join(input_path, "LULC_Products/"), "maxObs_clip.tif"),
		"pair_24_0" : os.path.join(os.path.join(input_path, "sea_pressure_24/"), "spressure_24_clip_%s.tif" %input_dates),
		"pair_inst_0" : os.path.join(os.path.join(input_path, "sea_pressure_inst/"), "spressure_i_clip_%s.tif" %input_dates),
		#"pair_inst" : "Pair_inst_%s.tif",
		"pre" : os.path.join(os.path.join(input_path,"precip_24/"), "precip_24_clip_%s.tif"),
		"hum_24" : os.path.join(os.path.join(input_path, "hum_24/"),"hum_24_clip_%s.tif" %input_dates),
		"hum_inst" : os.path.join(os.path.join(input_path, "hum_inst/"), "hum_i_clip_%s.tif" %input_dates),
		"tair_24" : os.path.join(os.path.join(input_path, "temp_24/"),"temp_24_clip_%s.tif" %input_dates),
		"tair_inst" : os.path.join(os.path.join(input_path, "temp_inst/"),"temp_i_clip_%s.tif" %input_dates),
		"tair_max_24" : os.path.join(os.path.join(input_path, "temp_24_max/"), "temp_max_clip_%s.tif" %input_dates),
		"tair_min_24" : os.path.join(os.path.join(input_path, "temp_24_min/"),"temp_min_clip_%s.tif" %input_dates),
		"tair_amp" : os.path.join(os.path.join(input_path, "temp_range"), "temp_range_clip.tif"),
		"wind_24" : os.path.join(os.path.join(input_path, "wind_24/"), "wind_clip%s.tif" %input_dates),
		"wind_inst" : os.path.join(os.path.join(input_path, "wind_inst/"), "wind_i_clip_%s.tif" %input_dates),
		"watCol_inst" : os.path.join(os.path.join(input_path, "watercol_inst/"), "watercol_i_clip_%s.tif" %input_dates),
		"trans_24" : os.path.join(os.path.join(input_path, "transmissivity_false/"), "transmissivity_clip.tif")
		}

		self.parDictIN = {
		"albedo" : os.path.join(os.path.join(input_path, "Albedo/"), "Albedo_%s.tif" %input_dates),
		"ndvi" : os.path.join(os.path.join(input_path, "Indecies/"), "Indices_%s.tif" %input_dates),
		"lst" : os.path.join(os.path.join(input_path, "LST/"), "LST_%s.tif" %input_dates),
		"time" : os.path.join(os.path.join(input_path, "Time/"), "solarTime.tif"),
		"lat" : os.path.join(os.path.join(input_path, "LatLon/"), "lat.tif"),
		"lon" : os.path.join(os.path.join(input_path, "LatLon/"), "lon.tif"),
		"dem" : os.path.join(os.path.join(input_path, "Z/"), "dem.tif"),
		"slope" : os.path.join(os.path.join(input_path, "Z/"), "slope.tif"),
		"aspect" : os.path.join(os.path.join(input_path, "Z/"), "aspect.tif"),
		"landMask" : os.path.join(os.path.join(input_path, "LandMask/"), "LandMask.tif"),
		"bulk" : os.path.join(os.path.join(input_path, "LULC_Products/"), "BSR.tif"),
		"maxObs" : os.path.join(os.path.join(input_path, "LULC_Products/"), "maxObs.tif"),
		"pair_24_0" : os.path.join(os.path.join(input_path, "sea_pressure_24/"), "spressure_24_%s.tif" %input_dates),
		"pair_inst_0" : os.path.join(os.path.join(input_path, "sea_pressure_inst/"), "spressure_i_%s.tif" %input_dates),
		#"pair_inst" : "Pair_inst_%s.tif",
		"pre" : os.path.join(os.path.join(input_path,"precip_24/"), "precip_24_%s.tif" %input_dates),
		"hum_24" : os.path.join(os.path.join(input_path, "hum_24/"),"hum_24_%s.tif" %input_dates),
		"hum_inst" : os.path.join(os.path.join(input_path, "hum_inst/"), "hum_i_%s.tif" %input_dates),
		"tair_24" : os.path.join(os.path.join(input_path, "temp_24/"),"temp_24_%s.tif" %input_dates),
		"tair_inst" : os.path.join(os.path.join(input_path, "temp_inst/"),"temp_i_%s.tif" %input_dates),
		"tair_max_24" : os.path.join(os.path.join(input_path, "temp_24_max/"), "temp_max_%s.tif" %input_dates),
		"tair_min_24" : os.path.join(os.path.join(input_path, "temp_24_min"), "temp_min_%s.tif" %input_dates),
		"tair_amp" : os.path.join(os.path.join(input_path, "temp_range"), "temp_range.tif"),
		"wind_24" : os.path.join(os.path.join(input_path, "wind_24/"), "wind_%s.tif" %input_dates),
		"wind_inst" : os.path.join(os.path.join(input_path, "wind_inst/"), "wind_i_%s.tif" %input_dates),
		"watCol_inst" : os.path.join(os.path.join(input_path, "watercol_inst/"), "watercol_i_%s.tif" %input_dates),
		"trans_24" : os.path.join(os.path.join(input_path, "transmissivity_false/"), "transmissivity.tif")
		}

		self.parDictOUT = {
		"vc" : os.path.join(output_path, "vc_%s.tif"),
		"lai" : os.path.join(output_path, "lai_%s.tif" %input_dates),
		"lai_eff" : os.path.join(output_path, "lai_eff_%s.tif" %input_dates),
		"sf_soil" : os.path.join(output_path, "sf_soil_%s.tif" %input_dates),
		"lat" : os.path.join(output_path, "lat_%s.tif" %input_dates),
		"slope" : os.path.join(output_path, "slope_%s.tif" %input_dates),
		"aspect" : os.path.join(output_path, "aspect_%s.tif" %input_dates),
		"ra_24_toa" : os.path.join(output_path, "ra_24_toa_%s.tif" %input_dates),
		"ws": os.path.join(output_path, "ws_%s.tif" %input_dates),
		"diffusion_index" : os.path.join(output_path, "diffusion_index_%s.tif" %input_dates),
		"ra_24" : os.path.join(output_path, "ra_24_%s.tif" %input_dates),
		"stress_rad" : os.path.join(output_path, "stress_rad_%s.tif" %input_dates),
		"p_air_24" : os.path.join(output_path, "_%s.tif" %input_dates),
		"vp_24" : os.path.join(output_path, "vp_24_%s.tif" %input_dates),
		"svp_24" : os.path.join(output_path, "svp_24_%s.tif" %input_dates),
		"vpd_24" : os.path.join(output_path, "vpd_24_%s.tif" %input_dates),
		"stress_vpd" : os.path.join(output_path, "stress_vpd_%s.tif" %input_dates),
		"stress_temp" : os.path.join(output_path, "stress_temp_%s.tif" %input_dates),
		"r_canopy_0" : os.path.join(output_path, "r_canopy_0_%s.tif" %input_dates),
		"t_air_k_24" : os.path.join(output_path, "t_air_k_24_%s.tif" %input_dates),
		"l_net" : os.path.join(output_path, "l_net_%s.tif" %input_dates),
		"int_mm" : os.path.join(output_path, "int_mm_%s.tif" %input_dates),
		"lh_24" : os.path.join(output_path, "lh_24_%s.tif" %input_dates),
		"int_wm2" : os.path.join(output_path, "int_wm2_%s.tif" %input_dates),
		"rn_24" : os.path.join(output_path, "rn_24_%s.tif" %input_dates),
		"rn_24_canopy" : os.path.join(output_path, "rn_24_canopy_%s.tif" %input_dates),
		"t_air_k_i" : os.path.join(output_path, "t_air_k_i_%s.tif" %input_dates),
		"vp_i" : os.path.join(output_path, "vp_i_%s.tif" %input_dates),
		"ad_moist_i" : os.path.join(output_path, "ad_moist_i_%s.tif" %input_dates),
		"ad_dry_i" : os.path.join(output_path, "ad_dry_i_%s.tif" %input_dates),
		"ad_i" : os.path.join(output_path, "ad_i_%s.tif" %input_dates),
		"u_b_i_bare" : os.path.join(output_path, "u_b_i_bare_%s.tif" %input_dates),
		"lon" : os.path.join(output_path, "lon_%s.tif" %input_dates),
		"ha" : os.path.join(output_path, "ha_%s.tif" %input_dates),
		"ied" : os.path.join(output_path, "ied_%s.tif" %input_dates),
		"h0" : os.path.join(output_path, "h0_%s.tif" %input_dates),
		"h0ref" : os.path.join(output_path, "h0ref_%s.tif" %input_dates),
		"m" : os.path.join(output_path, "m_%s.tif" %input_dates),
		"rotm" : os.path.join(output_path, "rotm_%s.tif" %input_dates),
		"Tl2"  : os.path.join(output_path, "Tl2_%s.tif" %input_dates),
		"B0c" : os.path.join(output_path, "B0c_%s.tif" %input_dates),
		"Bhc" : os.path.join(output_path, "Bhc_%s.tif" %input_dates),
		"Dhc" : os.path.join(output_path, "Dhc_%s.tif" %input_dates),
		"ra_hor_clear_i" : os.path.join(output_path, "ra_hor_clear_i_%s.tif" %input_dates),
		"emiss_atm_i"  : os.path.join(output_path, "emiss_atm_i_%s.tif" %input_dates),
		"rn_bare" : os.path.join(output_path, "rn_bare_%s.tif" %input_dates),
		"rn_full" : os.path.join(output_path, "rn_full_%s.tif" %input_dates),
		"u_b_i_full" : os.path.join(output_path, "u_b_i_full_%s.tif" %input_dates),
		"u_star_i_bare" : os.path.join(output_path, "u_star_i_bare_%s.tif" %input_dates),
		"u_star_i_full" : os.path.join(output_path, "u_star_i_full_%s.tif" %input_dates),
		"u_i_soil" : os.path.join(output_path, "u_i_soil_%s.tif" %input_dates),
		"ras"  : os.path.join(output_path, "ras_%s.tif" %input_dates),
		"raa" : os.path.join(output_path, "raa_%s.tif" %input_dates),
		"rac" : os.path.join(output_path, "rac_%s.tif" %input_dates),
		"t_max_bare" : os.path.join(output_path, "t_max_bare_%s.tif" %input_dates),
		"t_max_full" : os.path.join(output_path, "t_max_full_%s.tif" %input_dates),
		"w_i" : os.path.join(output_path, "w_i_%s.tif" %input_dates),
		"t_dew_i" : os.path.join(output_path, "t_dew_i_%s.tif" %input_dates),
		"t_wet_i" : os.path.join(output_path, "t_wet_i_%s.tif" %input_dates),
		"t_wet_k_i" : os.path.join(output_path, "t_wet_k_i_%s.tif" %input_dates),
		"lst_max"  : os.path.join(output_path, "lst_max_%s.tif" %input_dates),
		"se_root" : os.path.join(output_path, "se_root_%s.tif" %input_dates),
		"stress_moist" : os.path.join(output_path, "stress_moist_%s.tif" %input_dates),
		"r_canopy_0" : os.path.join(output_path, "r_canopy_0_%s.tif" %input_dates),
		"r_canopy" : os.path.join(output_path, "r_canopy_%s.tif" %input_dates),
		"z_obst" : os.path.join(output_path, "z_obst_%s.tif" %input_dates),
		"z_oro"  : os.path.join(output_path, "z_oro_%s.tif" %input_dates),
		"z0m" : os.path.join(output_path, "z0m_%s.tif" %input_dates),
		"ra_canopy_init" : os.path.join(output_path, "ra_canopy_init_%s.tif" %input_dates),
		"u_b_24" : os.path.join(output_path, "u_b_24_%s.tif" %input_dates),
		"disp" : os.path.join(output_path, "disp_%s.tif" %input_dates),
		"u_star_24_init" : os.path.join(output_path, "u_star_24_init_%s.tif" %input_dates),
		"ad_dry_24" : os.path.join(output_path, "ad_dry_24_%s.tif" %input_dates),
		"ad_moist_24" : os.path.join(output_path, "ad_moist_24_%s.tif" %input_dates),
		"ad_24" : os.path.join(output_path, "ad_24_%s.tif" %input_dates),
		"psy_24" : os.path.join(output_path, "psy_24_%s.tif" %input_dates),
		"ssvp_24" : os.path.join(output_path, "ssvp_24_%s.tif" %input_dates),
		"t_24_init" : os.path.join(output_path, "t_24_init_%s.tif" %input_dates),
		"h_canopy_24_init" : os.path.join(output_path, "h_canopy_24_init_%s.tif" %input_dates),
		"t_24" : os.path.join(output_path, "t_24_%s.tif" %input_dates),
		"t_24_mm" : os.path.join(output_path, "t_24_mm_%s.tif" %input_dates),
		"sf_soil" : os.path.join(output_path, "sf_soil_%s.tif" %input_dates),
		"rn_24_soil": os.path.join(output_path, "rn_24_soil_%s.tif" %input_dates),
		"r_soil" : os.path.join(output_path, "r_soil_%s.tif" %input_dates),
		"ra_soil_init" : os.path.join(output_path, "ra_soil_init_%s.tif" %input_dates),
		"u_b_24" : os.path.join(output_path, "u_b_24_%s.tif" %input_dates),
		"u_star_24_soil_init" : os.path.join(output_path, "u_star_24_soil_init_%s.tif" %input_dates),
		"g0_bs" : os.path.join(output_path, "g0_bs_%s.tif" %input_dates),
		"g0_24" : os.path.join(output_path, "g0_24_%s.tif" %input_dates),
		"e_24_init" : os.path.join(output_path, "e_24_init_%s.tif" %input_dates),
		"h_soil_24_init" : os.path.join(output_path, "h_soil_24_init_%s.tif" %input_dates),
		"e_24" : os.path.join(output_path, "e_24_%s.tif" %input_dates),
		"e_24_mm" : os.path.join(output_path, "e_24_mm_%s.tif" %input_dates),
		"et_24_mm" : os.path.join(output_path, "et_24_mm_%s.tif" %input_dates),
		"rn_24_grass" : os.path.join(output_path, "rn_24_grass_%s.tif" %input_dates),
		"et_ref_24" : os.path.join(output_path, "et_ref_24_%s.tif" %input_dates),
		"et_ref_24_mm" : os.path.join(output_path, "et_ref_24_mm_%s.tif" %input_dates)
		}

	def getFilePathIN(self, key):
		return self.parDictIN.get(key)
	
	def getKeysIN(self):
		return self.parDictIN.keys()

	def getFilePathOUT(self, key):
		return self.parDictOUT.get(key)
	
	def getClipPathIN(self, key):
		return self.parDictClip.get(key)

	def getExtent(self):
		return self.extent

