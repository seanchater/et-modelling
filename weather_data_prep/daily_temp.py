import pandas
import numpy

# read in full weather data csv
filename = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data\ws\Hourly_Data\Climate_Data.csv"
df = pandas.read_csv(filename)

# check
print(df.describe())
print(df.head())

# select only the data that we are processing (temperature for this version)
temp_df = df[["id", "date", "station_id", "time", "temperature_mean", "latitude", "longitude"]]

print(temp_df.head())

# we cant use null values so remove (cant replace with zero becuase the temp might not actually be zero)
temp = temp_df[temp_df["temperature_mean"].notna()]

# create a new col containing month that will be used to group stuff 
temp["month"] = temp.apply(lambda row: row.date.split("-")[1], axis=1)
temp.head()

# first select all observations from a single month (eg: 01)
d = {
	1: "01",
	2: "02",
	3: "03", 
	4: "04",
	5: "05", 
	6: "06",
	7: "07",
	8: "08",
	9: "09",
	10: "10",
	11: "11",
	12: "12"
}

for i in range(1,13):

	month_obs = temp[temp["month"] == d[i]]

	# only use temperatures within an acceptable range [-10, 50]
	temp_valid = month_obs[(month_obs["temperature_mean"] > -10) & (month_obs["temperature_mean"] < 50)]

	# then group by station (position: lat/lon) to get average temp of station for month
	temp_group = temp_valid.groupby(["latitude", "longitude"])["temperature_mean"].mean()

	# export to csv
	outname = "{}_av_temp.csv".format(d[i])
	temp_group.to_csv(outname)