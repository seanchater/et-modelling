# TODO: import data
'''
Data needed:

INPUT						TYPE			SENSOR 			PRODUCT 			COMMENT

Precipitation 				Model 							CHIRPS v2
															CHIRP


Surface albedo 				Sensor 			MODIS 			MOD09GA,
											MOD09GQ


Weather data (temp,
specific humidity, wind
speed, air pressure,
aerosol optical depth)		Model							MERRA/GEOS5			MERRA used prior to start of GEOS-5 (21-2-2014)


NDVI 						Sensor 			MODIS			MOD09GQ


Land Surface
Temperature					Sensor 			MODIS 			MOD11A1,

															MYD11A1				Used to derive Soil Moisture Stress

Elevation, slope and
aspect						Static 							SRTM 				Elevation, slope and aspect are derived from the DEM


Transmissivity 				Model 			MSG 								Transmissivity is derived from MSG shortwave radiation products


Land Cover 					Static 			WaPOR L1
											Land Cover 
											Classification						If L1 LC was not yet available, preliminary LC obtained from ESA GlobCover
'''

# TODO: preproces and prepare inputs
'''
If using ASMRE data then might need to do resampling coz 25km res
'''

# TODO: calculate E & T using PM equations

# TODO: produce output maps (raster grids probably)