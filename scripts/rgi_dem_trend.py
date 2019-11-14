#! /usr/bin/env python

"""
Script to create DEMStack objects for list of DEMs over RGI polygons

Run dem_align_post.py and rgi_dem_trend_prep.sh first

To do:
Should run in parallel with multiprocessing, similar to mb_parallel

"""

import os
import sys
import re
from osgeo import ogr
from pygeotools.lib import geolib, warplib, malib

#Add proper argument parsing

dem_index_fn = sys.argv[1]

#Minimum glacier area (km2)
#min_glac_area = 0.1 
min_glac_area = float(sys.argv[2])
max_glac_area = float(sys.argv[3])

#Buffer this distance (m) around RGI polygon for stats
buffer_m = 1000

#Use 2 threads for the polygons with smaller area
if min_glac_area >= 2.0: 
    n_threads = 14 
else:
    n_threads = 2 

#Most DEMs are 32 m
#res='max'
res=30
#res=8

#Minimum number of samples
min_dem_count = 5
#min_dem_count = 3
#Min to max timestamp difference (days)
#5 years
min_dt_ptp = 1825
#min_dt_ptp = 1095

topdir='/nobackup/deshean/'

#demdir = os.path.join(topdir, 'hma/aster/dsm')
#demdir = os.path.join(topdir, 'hma/dem_coreg')
demdir = os.path.join(topdir, 'hma/combined_aster_wv')
#demdir = os.path.join(topdir, 'hma/aster/dsm/dem_align_ASTERonly')

os.chdir(demdir)
stackdir = os.path.splitext(dem_index_fn)[0]+'_stack'
if not os.path.exists(stackdir):
    os.makedirs(stackdir)

#Output log for each file, avoid filling PBS spool
#https://www.nas.nasa.gov/hecc/support/kb/avoiding-job-failure-from-overfilling-pbsspool_183.html
logdir=os.path.join(stackdir, 'log')
if not os.path.exists(logdir):
    os.makedirs(logdir)

def rgi_name(feat, glacname_fieldname='Name', glacnum_fieldname='RGIId'):
    #glacname = feat.GetField(feat.GetFieldIndex(glacname_fieldname))
    glacname = feat.GetField(glacname_fieldname)
    if glacname is None:
        glacname = ""
    else:
        #RGI has some nonstandard characters
        #This worked in Python2
        #glacname = glacname.decode('unicode_escape').encode('ascii','ignore')
        #This should be universal
        #glacname = re.sub(r'[^\x00-\x7f]',r'', glacname)
        glacname = re.sub(r'\W+', '', glacname)
        glacname = glacname.replace(" ", "")
        glacname = glacname.replace("_", "")
        glacname = glacname.replace("/", "")

    #glacnum = feat.GetField(feat.GetFieldIndex(glacnum_fieldname))
    glacnum = feat.GetField(glacnum_fieldname)
    fn = feat.GetDefnRef().GetName()
    #RGIId (String) = RGI50-01.00004
    #glacnum = '%0.5f' % float(glacnum.split('-')[-1])
    glacnum = float(glacnum.split('-')[-1])

    if glacname:
        feat_fn = "%08.5f_%s" % (glacnum, glacname)
    else:
        feat_fn = "%08.5f" % glacnum
    return feat_fn

glac_shp_fn = os.path.join(topdir,'data/rgi60/regions/rgi60_merge_HMA_aea.shp')
glac_shp_ds = ogr.Open(glac_shp_fn, 0)
glac_shp_lyr = glac_shp_ds.GetLayer()
#This should be contained in features
glac_shp_srs = glac_shp_lyr.GetSpatialRef()
feat_count = glac_shp_lyr.GetFeatureCount()
print("Input glacier polygon count: %i" % feat_count)

#Area filter
glac_shp_lyr.SetAttributeFilter("Area > %s and Area < %s" % (min_glac_area, max_glac_area))
feat_count = glac_shp_lyr.GetFeatureCount()
print("Min. Area filter glacier polygon count: %i" % feat_count)
glac_shp_lyr.ResetReading()

dem_index_ds = ogr.Open(dem_index_fn, 0)
dem_index_lyr = dem_index_ds.GetLayer()
#This should be contained in features
dem_index_srs = dem_index_lyr.GetSpatialRef()
feat_count = dem_index_lyr.GetFeatureCount()
print("Input dem count: %i" % feat_count)

#cmd_fn = 'dem_stack_cmd.sh'
cmd_fn = os.path.splitext(dem_index_fn)[0]+'_%s-%s_km2_stack_cmd.sh' % (min_glac_area, max_glac_area)
f = open(cmd_fn, "w") 

glac_dict = None
#glac_dict = ['15.03473', '15.03733', '15.10070', '15.09991']

for n, feat in enumerate(glac_shp_lyr):
    #glac_geom_orig = geolib.geom_dup(feat.GetGeometryRef())
    feat_fn = rgi_name(feat)
    if glac_dict is not None:
        if not feat_fn in glac_dict:
            continue
    print(n, feat_fn)
    glac_geom = geolib.geom_dup(feat.GetGeometryRef())
    glac_geom_extent = geolib.geom_extent(glac_geom)
    #print(glac_geom_extent)
    #Should buffer by ~1 km here, preserve surrounding pixels for uncertainty analysis
    glac_geom = glac_geom.Buffer(buffer_m)
    glac_geom_extent = geolib.geom_extent(glac_geom)
    #print(glac_geom_extent)
    #glac_geom_extent = geolib.pad_extent(glac_geom_extent, width=1000)

    #Spatial filter
    dem_index_lyr.SetSpatialFilter(glac_geom)
    dem_count = dem_index_lyr.GetFeatureCount()
    print("dem count after spatial filter: %i" % dem_count)
    if dem_count > min_dem_count:
        fn_list = []
        for dem_feat in dem_index_lyr:
            #Only 1 field from gdaltindex, 'location'
            #Can have issues with commands that are too long with full path
            #fn = os.path.join(demdir, dem_feat.GetField(0))
            #This is relative path, must be in correct directory when running make_stack.py
            fn = dem_feat.GetField(0)
            fn_list.append(fn)
        fn_list.sort() 
        #Hack to deal with long filenames
        outdir=os.path.join(stackdir, feat_fn)
        #stack_fn='%s_%s_%s_stack.npz' % (feat_fn[0:8], os.path.split(fn_list[0])[-1].split('_DEM_')[0], os.path.split(fn_list[-1])[-1].split('_DEM_')[0])
        stack_fn='%s_%s.npz' % (feat_fn[0:8], os.path.split(stackdir)[-1])

        logfile = os.path.join(logdir, os.path.splitext(stack_fn)[0]+'.log')

        #Create file with commands to make stacks
        #Run this output file with GNU parallel `parallel < dem_stack_cmd.sh`
        #For glaciers with smaller areas, use -j 28; For glaciers with larger areas, use -j 4
        cmd='make_stack.py --med --trend --robust'
        cmd+=' -n_cpu %i -outdir %s -stack_fn %s -tr %s -te "%s" -t_srs "%s" -min_n %i -min_dt_ptp %f %s > %s 2>&1 \n' % \
                (n_threads, outdir, os.path.join(outdir, stack_fn), res, ' '.join(str(i) for i in glac_geom_extent), \
                dem_index_srs.ExportToProj4(), min_dem_count, min_dt_ptp, ' '.join(fn_list), logfile)
        f.write(cmd)

        #Generate stacks serially
        #stack = malib.DEMStack(fn_list, outdir=os.path.join(stackdir, feat_fn), \
        #        res='max', extent=glac_geom_extent, srs=dem_index_srs, mask_geom=glac_geom, \
        #        trend=True, robust=True, n_thresh=min_dem_count, min_dt_ptp=min_dt_ptp)

        #if stack.ma_stack is not None:
            #sys.exit()
            #glac_geom_mask = geolib.geom2mask(glac_geom, stack.get_ds())
            #ds_list = warplib.memwarp_multi_fn(fn_list, res='max', extent=glac_geom_extent)

    dem_index_lyr.ResetReading()

f.close()

#Generate plots
#hs.sh */*med.tif
#for i in */*trend.tif; do  imviewer.py -clim -5 5 -cmap RdYlBu -label 'Trend (m/yr)' -of png -overlay $(echo $i | sed 's/_trend/_med_hs_az315/') $i -scalebar -outsize 8 8 -dpi 100; done
#for i in */*count.tif; do imviewer.py -clim 0 20 -cmap inferno -label 'Count' -of png -overlay $(echo $i | sed 's/_count/_med_hs_az315/') $i -scalebar -outsize 8 8 -dpi 100; done
#for i in */*[0-9]_nmad.tif; do imviewer.py -clim 0 30 -label 'Std (m)' -of png -overlay $(echo $i | sed 's/_std/_med_hs_az315/') $i -scalebar -outsize 8 8 -dpi 100; done
