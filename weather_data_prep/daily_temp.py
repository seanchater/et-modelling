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
obs_jan = temp[temp["month"] == "01"]

# only use temperatures within an acceptable range [-10, 50]
temp = temp[(temp["temperature_mean"] > -10) & (temp["temperature_mean"] < 50)]

# then group by station (position: lat/lon) to get average temp of station for month
temp_group = temp.groupby(["latitude", "longitude"])["temperature_mean"].mean()

# export to csv
temp_group.to_csv("jan_av_temp.csv")