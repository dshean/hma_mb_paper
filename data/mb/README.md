# HMA mass balance data

These csv files contain a row of floating point values derived for each RGI v6.0 glacier polygon. 

Please see script used to generate for more detail: https://github.com/dshean/hma_mb_paper/blob/master/scripts/mb_parallel.py

Variable names written to csv: https://github.com/dshean/hma_mb_paper/blob/master/scripts/mb_parallel.py#L1247

## Columns
The first line of the csv contains a header:  
`RGIId,x,y,z_med,z_min,z_max,z_slope,z_aspect,dhdt_ma,dhdt_ma_sigma,mb_mwea,mb_mwea_sigma,area_m2,mb_m3wea,mb_m3wea_sigma,t1,t2,dt,valid_area_perc,H_m,debris_m,perc_debris,perc_pond,perc_clean,vm_ma`

Sample row (Ngozumpa Glacier, Nepal):  
`15.03474,159420.734,-893584.333,5492.891,4903.261,7219.088,13.396,181.068,-0.552,0.232,-0.469,0.201,12919247.330,-6064497.893,2653695.223,2000.412,2018.412,18.000,99.968,94.255,0.281,42.958,0.719,56.323,6.321`

### Standard parameters
|Parameter|Description|Units|Example|
|---|---|---|---|
|RGIId|RGI glacier polygon ID|15.03474| 
|x|glacier polygon centroid x-coordinate in custom Albers-Equal-Area coordinate system. See text for details, PROJ4 string `'+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '`)|meters|159420.734| 
|y|glacier polygon centroid y-coordinate in custom Albers-Equal-Area coordinate system. See text for details, PROJ4 string `'+proj=aea +lat_1=36 +lat_2=49 +lat_0=43 +lon_0=-115 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs '`)|meters|-893584.333| 
|z_med|median glacier elevation from ASTER/WV DEM composite|meters above WGS84 ellipsoid|5492.891| 
|z_min|minimum glacier elevation from ASTER/WV DEM composite|meters above WGS84 ellipsoid|4903.261| 
|z_max|maximum glacier elevation from ASTER/WV DEM composite|meters above WGS84 ellipsoid|7219.088| 
|z_slope|median glacier slope from ASTER/WV DEM composite|degrees|13.396| 
|z_aspect|median glacier aspect from ASTER/WV DEM composite|degrees (clockwise from North)|181.068| 
|dhdt_ma|mean glacier elevation change (see text for details)|meters per year|-0.552| 
|dhdt_ma_sigma|glacier elevation change uncertainty (see text for details)|meters per year|0.232| 
|mb_mwea|mean specific glacier mass balance (see text for details)|meters water equivalent per year|-0.469| 
|mb_mwea_sigma|specific glacier mass balance uncertainty (see text for details)|meters water equivalent per year|0.201| 
|area_m2|glacier polygon valid dh/dt area|m^2|12919247.330| 
|mb_m3wea|total glacier mass balance (see text for details)|cubic meters water equivalent per year|-6064497.893|
|mb_m3wea_sigma|total glacier mass balance uncertainty (see text for details)|cubic meters water equivalent per year|2653695.223|
|t1|Begin date|decimal year|2000.412|
|t2|End date|decimal year|2018.412|
|dt|Time difference|decimal year|18.000|
|valid_area_perc|glacier polygon valid dh/dt area as percentage of total polygon area|%|99.968| 

### Auxiliary parameters
Included for glaciers with grids available from external sources, glaciers larger than 2 km^2

* Ice thickness: Updated to RGI v6.0 by Matthias Huss, downloaded March 2018, based on Huss and Farinotti (2012)
* Debris classification: Kraaijenbrink et al. (2017), http://mountainhydrology.org/data-nature-2017/
* Surface Velocity: Preliminary ITS_LIVE data, downloaded August 2018 (HMA_G0240_0000), https://its-live.jpl.nasa.gov/ 

|Parameter|Description|Units|Example|
|---|---|---|---|
|H_m|Mean ice thickness from Matthias Huss (2018)|meters|94.255|
|debris_m|Mean debris thickness from Kraaijenbrink et al. (2017)|meters|0.281|
|perc_debris|Percent of glacier polygon classified as debris, from Kraaijenbrink et al. (2017)|%|42.958|
|perc_pond|Percent of glacier polygon classified as pond, from Kraaijenbrink et al. (2017)|%|0.719|
|perc_clean|Percent of glacier polygon classified as clean ice, from Kraaijenbrink et al. (2017)|%|56.323|
|vm_ma|Mean glacier velocity magnitude from ITS_LIVE, Gardner et al. (2019)|meters per year|6.321|
