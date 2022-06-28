from typing import Dict


class PARAMS:
    def __init__(self, input_path, input_dates:None):

        import os

        self.parDict = {
            "albedo" : os.path.join(os.path.join(input_path, "Albedo_outputs/"), "Albedo_%s.tif" %input_dates),
            "ndvi" : os.path.join(os.path.join(input_path, "Indecies/"), "Indeces_%s.tif" %input_dates),
            "lst" : os.path.join(os.path.join(input_path, "LST/"), "LST_%s.tif" %input_dates),
            #"time" : os.path.join(input_path, "Time_%s.tif"),
            "lat" : os.path.join(os.path.join(input_path, "Lat_Long_Rasters/"), "Lat.tif"),
            "lon" : os.path.join(os.path.join(input_path, "Lat_Long_Rasters/"), "Lon.tif"),
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
            #"wind_24" : "wind_24_%s.tif",
            #"wind_inst" : "wind_inst_%s.tif",
            #"watCol_inst" : "wv_inst_%s.tif",
            "trans_24" : os.path.join(os.path.join(input_path, "transmissivity_false/"), "transmissivity.tif")
        }

    def getFileName(self, key):
        return self.parDict.get(key)
    
