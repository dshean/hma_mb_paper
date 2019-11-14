#! /usr/bin/env python

#Filter and interpolate stack dhdt to specified timestamps

import os
import sys
import numpy as np
from datetime import datetime
from pygeotools.lib import malib, timelib, iolib, filtlib

#Need to do proper argument parsing
#Accept filter size, dz diff, dt_list, etc

#These are from TanDEM-X masked
#Not sure why negative outliers
#zlim = (-1013,8388)
#Ellispoidal heights at sea level should be above -100 m
zlim = (-200,8388)
zlim_pad = 150

filter=True

#Window size for median and gauss filters (px)
size=3

#Input stack npz
stack_fn=sys.argv[1]

#User-provided timestamps with format YYYYMMDD
if len(sys.argv) > 2:
    dt_list_str = sys.argv[2:]
    if len(dt_list_str) == 1:
        dt_list_str = dt_list_str[0].split(' ')
    dt_list = [timelib.fn_getdatetime_list(dt_str)[0] for dt_str in dt_list_str]
else:
    #SRTM, then systematic ASTER timestamps
    dt_list = [datetime(2000,2,11), datetime(2000,5,31), datetime(2009,5,31), datetime(2018,5,31)]

#Use tif on disk if available
out_fn=os.path.splitext(stack_fn)[0]
#Otherwise load stack and compute trend/intercept if necessary

trend_fn=out_fn+'_trend.tif'
trend_ds = iolib.fn_getds(trend_fn)
trend = iolib.ds_getma(trend_ds)/365.25

intercept_fn=out_fn+'_intercept.tif'
#Hmmm, no 365.25 factor here.  Clean up here and in stack generation
intercept = iolib.fn_getma(intercept_fn)

#Determine local elevation min/max - use to throw out bad trends
med_fn = out_fn+'_med.tif'
if not os.path.exists(med_fn):
    med_fn = out_fn+'_mean.tif'
med = iolib.fn_getma(med_fn)
zlim = malib.calcperc(med, (0.01, 99.99))
zlim = (zlim[0] - zlim_pad, zlim[1] + zlim_pad)
med = None

#Can vectorize
#dt_list_o = timelib.dt2o(dt_list)
#z_list = trend*dt_list_o[:,None,None]+intercept

if filter:
    #Could remove outliers in trend at this phase
    print("Input pixel count: %s" % trend.count())

    #Outlier filter
    #This can remove valid pixels for larger glaciers (e.g. Baltoro) with negative dh/dt - could scale based on pixel count
    #trend_filt = filtlib.mad_fltr(trend, n=4)
    #trend_filt = filtlib.sigma_fltr(trend, n=3)
    #print("Output pixel count: %s" % trend_filt.count())

    #Remove islands
    #trend_filt = filtlib.remove_islands(trend, iterations=1)
    #Erode edges near nodata 
    trend_filt = filtlib.erode_edge(trend, iterations=1)
    print("Output pixel count: %s" % trend_filt.count())

    #Rolling median filter (remove noise) - can use a slightly larger window here
    trend_filt = filtlib.rolling_fltr(trend_filt, size=size, circular=True, origmask=True)
    print("Output pixel count: %s" % trend_filt.count())

    #Gaussian filter (smooth)
    #trend_filt = filtlib.gauss_fltr_astropy(trend_filt, size=size, origmask=True, fill_interior=True)
    trend_filt = filtlib.gauss_fltr_astropy(trend_filt, size=size)
    print("Output pixel count: %s" % trend_filt.count())

    trend_fn=out_fn+'_trend_%spx_filt.tif' % size
    print("Writing out: %s" % trend_fn)
    iolib.writeGTiff(trend_filt*365.25, trend_fn, trend_ds)

    #Update intercept using new filtered slope values
    #Need to update for different periods?
    #dt_pivot = timelib.mean_date(datetime(2000,5,31), datetime(2009,5,31))
    #dt_pivot = timelib.mean_date(datetime(2009,5,31), datetime(2018,5,31))
    dt_pivot = timelib.dt2o(datetime(2009, 5, 31))
    intercept_filt = dt_pivot * (trend - trend_filt) + intercept
    intercept_fn=out_fn+'_intercept_%spx_filt.tif' % size
    print("Writing out: %s" % intercept_fn)
    iolib.writeGTiff(intercept_filt*365.25, intercept_fn, trend_ds)

    trend = trend_filt
    intercept = intercept_filt

for dt in dt_list:
    dt_o = timelib.dt2o(dt)
    z = trend*dt_o+intercept
    #Remove any values outsize global limits
    #Could also do local range filter here
    z = filtlib.range_fltr(z, zlim)
    print("Output pixel count: %s" % z.count())
    interp_fn=os.path.splitext(trend_fn)[0]+'_%s.tif' % dt.strftime('%Y%m%d')
    print("Writing out: %s" % interp_fn)
    iolib.writeGTiff(z, interp_fn, trend_ds)
