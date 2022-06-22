import ee
ee.Authenticate()
ee.Initialize()


#Code starts here
#variable selection for LST and DEM
lst=ee.ImageCollection('MODIS/006/MOD11A1')

srtdem=ee.Image("CGIAR/SRTM90_V4")
elevation=srtdem.select('elevaion')
slope=ee.Terrain.slope(elevation)

#Date and filtering data
#start date is inclusive end date is exlusive 
i_date = '2021-01-01'
f_date = '2021-01-02'

lst=lst.select('LST_Day_1km').filterDate(i-date, f_date)
