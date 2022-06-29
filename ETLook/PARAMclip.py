class PARAMclip:
    def __init__(self, input_path, output_path, input_dates:None) -> None:
        
        import os

        self.extent = os.path.join(os.path.join(input_path, "3318DD_extent/"), "3318DD.shp")

        self.parDictIN = {
            "albedo" : os.path.join(os.path.join(input_path, "Albedo_outputs/"), "Albedo_%s.tif" %input_dates),
            "ndvi" : os.path.join(os.path.join(input_path, "Indecies/"), "Indeces_%s.tif" %input_dates),
            "lst" : os.path.join(os.path.join(input_path, "LST/"), "LST_%s.tif" %input_dates),
            "time" : os.path.join(os.path.join(input_path, "Time/"), "solarTime.tif"),
            "lat" : os.path.join(os.path.join(input_path, "Lat_Long_Rasters/"), "lat.tif"),
            "lon" : os.path.join(os.path.join(input_path, "Lat_Long_Rasters/"), "lon.tif"),
            "dem" : os.path.join(os.path.join(input_path, "DEM_Derviatives/"), "dem.tif"),
            "slope" : os.path.join(os.path.join(input_path, "DEM_Derviatives/"), "slope.tif"),
            "aspect" : os.path.join(os.path.join(input_path, "DEM_Derviatives/"), "aspect.tif"),
            "landMask" : os.path.join(os.path.join(input_path, "Land_Mask/"), "LandMask.tif"),
            "bulk" : os.path.join(os.path.join(input_path, "LULC_Products/"), "BSR_AOI.tif"),
            "maxObs" : os.path.join(os.path.join(input_path, "LULC_Products/"), "Max_height_AOI.tif"),
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
            "wind_24" : os.path.join(os.path.join(input_path, "Wind_IDW/"), "Wind_IDW_24.tif"),
            #"wind_inst" : "wind_inst_%s.tif",
            #"watCol_inst" : "wv_inst_%s.tif",
            "trans_24" : os.path.join(os.path.join(input_path, "transmissivity_false/"), "transmissivity.tif")
        }

        self.parDictClip = {
            "albedo" : os.path.join(os.path.join(input_path, "Albedo_outputs/"), "albedo_clip_%s.tif" %input_dates),
            "ndvi" : os.path.join(os.path.join(input_path, "Indecies/"), "ndvi_clip_%s.tif" %input_dates),
            "lst" : os.path.join(os.path.join(input_path, "LST/"), "lst_clip_%s.tif" %input_dates),
            "time" : os.path.join(os.path.join(input_path, "Time/"), "solarTime_clip.tif"),
            "lat" : os.path.join(os.path.join(input_path, "LatLong/"), "lat_clip.tif"),
            "lon" : os.path.join(os.path.join(input_path, "LatLong/"), "lon_clip.tif"),
            "dem" : os.path.join(os.path.join(input_path, "DEM_/"), "dem_clip.tif"),
            "slope" : os.path.join(os.path.join(input_path, "DEM_/"), "slope_clip.tif"),
            "aspect" : os.path.join(os.path.join(input_path, "DEM_/"), "aspect_clip.tif"),
            "landMask" : os.path.join(os.path.join(input_path, "LandMask/"), "landMask_clip.tif"),
            "bulk" : os.path.join(os.path.join(input_path, "LULC/"), "bulkSt_clip.tif"),
            "maxObs" : os.path.join(os.path.join(input_path, "LULC/"), "maxObs_clip.tif"),
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
            "wind_24" : os.path.join(os.path.join(input_path, "Wind/"), "wind_24_clip.tif"),
            #"wind_inst" : "wind_inst_%s.tif",
            #"watCol_inst" : "wv_inst_%s.tif",
            "trans_24" : os.path.join(os.path.join(input_path, "transmissivity/"), "transmissivity_clip.tif")
        }

    def getFilePathIN(self, key):
        return self.parDictIN.get(key)
    
    def getFilePathOUT(self, key):
        return self.parDictOUT.get(key)

    def getExtent(self):
        return self.extent


