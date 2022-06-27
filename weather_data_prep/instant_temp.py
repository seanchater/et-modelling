import pandas

filename = r"C:\Users\seanc\Documents\SU\2022_hons\716\et\etlook\input_data\ws\Hourly_Data\Climate_Data.csv"
df = pandas.read_csv(filename)

print(df.describe())

temp_df = df[["id", "date", "station_id", "time", "temperature_mean", "latitude", "longitude"]]

# print(temp_df.head())
# print(temp_df.shape)
temp = temp_df.fillna(0)
temp = temp_df[temp_df["temperature_mean"].notna()]

temp["month"] = temp.apply(lambda row: row.date.split("-")[1], axis=1)
temp["day"] = temp.apply(lambda row: row.date.split("-")[2], axis=1)

# select a single time to get temperatures (eg: 600 = 06:00am)
time_obs = temp[temp["time"] == 1000]
# print(time_obs.shape)
# print(time_obs.head(10))

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
		
	month_obs = time_obs[time_obs["month"] == d[i]]

	# only use temperatures within an acceptable range [-10, 50]
	temp_valid = month_obs[(month_obs["temperature_mean"] > -10) & (month_obs["temperature_mean"] < 50)]

	# get median of each station for this current month
	temp_group = temp_valid.groupby(["latitude", "longitude"])["temperature_mean"].median()
	# print(temp_group.shape)
	# print(temp_group.head())

	# export to csv
	outname = "{}_inst_temp.csv".format(d[i])
	temp_group.to_csv(outname)