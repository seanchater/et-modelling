# choose one of the two options below
#ra_24 = ETLook.solar_radiation.daily_solar_radiation_flat(ra_24_toa_flat, trans_24)
ra_24 = ETLook.solar_radiation.daily_total_solar_radiation(ra_24_toa, ra_24_toa_flat, diffusion_index, trans_24)
stress_rad = ETLook.stress.stress_radiation(ra_24)
p_air_24 = ETLook.meteo.air_pressure_daily(z, p_air_0_24)
vp_24 = ETLook.meteo.vapour_pressure_from_specific_humidity_daily(qv_24, p_air_24)
svp_24 = ETLook.meteo.saturated_vapour_pressure_average(
            ETLook.meteo.saturated_vapour_pressure_maximum(t_air_max_24),
            ETLook.meteo.saturated_vapour_pressure_minimum(t_air_min_24))
vpd_24 = ETLook.meteo.vapour_pressure_deficit_daily(svp_24, vp_24)
stress_vpd = ETLook.stress.stress_vpd(vpd_24, vpd_slope)
stress_temp = ETLook.stress.stress_temperature(t_air_24, t_opt, t_min, t_max)
r_canopy_0 = ETLook.resistance.atmospheric_canopy_resistance(lai_eff, stress_rad, stress_vpd, stress_temp, rs_min, rcan_max)