# et-modelling
ET Look model automation project for GIT 716.

- [et-modelling](#et-modelling)
  - [Overview of Input and Output Data](#overview-of-input-and-output-data)
    - [Input Files:](#input-files)
    - [Output Files:](#output-files)
    - [File Structure:](#file-structure)
  - [Conceptual Model of Code:](#conceptual-model-of-code)

## Overview of Input and Output Data

### Input Files:

Inputs | Range Required
:---| :---
Albedo|Yes
NDVI|Yes
LST (Land Surface Temperature)|Yes
Time (Solar Time of Capture)| No
Lat (Grid of Latitude values)|No
Lon (Grid of Longitude values)|No
DEM|No
Slope|No
Aspect|No
Land Mask|No
Bulk Stomatal|No
Max Obstacle Height|No
Air Pressure (Daily)|Yes
Air Pressure (Instantaneous)|Yes
Precipitation|No
Humidity (Daily)|Yes
Humidity (Instantaneous)|Yes
Air Temperature (Daily)|Yes
Air Temperature (Instantaneous)|Yes
Air Temperature (Daily max)|Yes
Air Temperature (Daily min)|Yes
Air Temperature (Amplitude)|No
Wind (Daily)|Yes
Wind (Instantaneous)|Yes
Transpiration|No

### Output Files:

Outputs | Saved By Default
:---| :---
vc |No
lai|No
lai_eff|No
sf_soil|No
lat|No
slope|No
aspect|No
ra_24_toa|No
ws|No
diffusion_index|No
ra_24|No
stress_rad|No
p_air_24|No
vp_24|No
svp_24|No
vpd_24|No
stress_vpd|No
stress_temp|No
r_canopy_0|No
t_air_k_24|No
l_net|No
int_mm|No
lh_24|No
int_wm2|No
rn_24|No
rn_24_canopy|No
t_air_k_i|No
vp_i|No
ad_moist_i|No
ad_dry_i|No
ad_i|No
u_b_i_bare|No
lon|No
ha|No
ied|No
h0|No
h0ref|No
m|No
rotm|No
Tl2|No
B0c|No
Bhc|No
Dhc|No
ra_hor_clear_i|No
emiss_atm_i|No
rn_bare|No
rn_full|No
u_b_i_full|No
u_star_i_bare|No
u_star_i_full|No
u_i_soil|No
ras|No
raa|No
rac|No
t_max_bare|No
t_max_full|No
w_i|No
t_dew_i|No
t_wet_i|No
t_wet_k_i|No
lst_max|No
se_root|No
stress_moist|No
r_canopy_0|No
r_canopy|No
z_obst|No
z_oro|No
z0m|No
ra_canopy_init|No
u_b_24|No
disp|No
u_star_24_init|No
ad_dry_24|No
ad_moist_24|No
ad_24|No
psy_24|No
ssvp_24|No
t_24_init|No
h_canopy_24_init|No
t_24|No
t_24_mm|No
sf_soil|No
rn_24_soil|No
r_soil|No
ra_soil_init|No
u_b_24|No
u_star_24_soil_init|No
g0_bs|No
g0_24|No
e_24_init|No
h_soil_24_init|No
e_24|No
e_24_mm|Yes
et_24_mm|Yes
rn_24_grass|No
et_ref_24|Yes
et_ref_24_mm|Yes

### File Structure:

Default Path Stucture| | | |
:---|:---|:---|:---
parent| | | |
 *|Data|||
*|~|input_Data| |
*|~|output| |
*|et-modelling||| 
*|~|ETLook||
*||~|Date-Range CSV|  

*By default, the code will follow the following file structure:*

## Conceptual Model of Code:

A conceptual model showcasing a flow chart of how all the algorithms integrate can be seen [here.](https://lucid.app/documents/view/a9214295-baff-4df2-bb31-5204af40405c)

Please see the flow chart for understanding the process.

It should also be noted that a lot of or code was referenced from the original [pyWAPOR](https://github.com/DHI-GRAS/wapor-et-look) code, which can also be found on GitHub.