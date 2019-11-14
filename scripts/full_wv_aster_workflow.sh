#! /bin/bash

#High-level wrapper and notes for DEM processing, geodetic mass balance and mosaic generation
#Feb 2019 

aster=/nobackup/deshean/hma/aster/dsm
wv=/nobackup/deshean/hma/dem_coreg/*track

#Filtering
#Abs diff +/-200 m to SRTM void-filled

#Prepare reference DEM
#~/src/tandemx/tandemx_eval.ipynb
~/src/tandemx/tandemx_proc.sh

#Co-registration
#refdem=/nobackup/deshean/data/tandemx/hma/mos_warp_erode/TDM1_DEM_90m_hma_DEM_masked_aea.vrt
#This runs dem_align.py on each DEM
#Run separately for ASTER and WV/GE
qsub ~/src/demcoreg/demcoreg/dem_align_parallel.pbs

#Review results, filter outliers, generate figures, remove bad products, generate combined?
#Start jupyterhub on node
~/src/demcoreg/demcoreg/dem_align_post.ipynb
#Old filtering/plotting script: ~/src/gmbtools/gmbtools/dem_align_post.py

#Generate indices and commands for each dataset for different time periods (gdaltindex, ogr_merge, rgi_dem_trend.py)
#Creates combined directory 'combined_aster_wv', years subdirectory with symlinks to align.tif products
~/src/gmbtools/gmbtools/rgi_dem_trend_prep.sh

#Generate stacks, produces trend.tif products
#Review/edit dem_post_parallel.pbs
qsub ~/src/gmbtools/gmbtools/dem_post_parallel.pbs

#Interpolate to end dates for each stack using stack_interp.py, generates vrt and retiles
#stack_interp also filters outliers, cleans up trend maps 
#Review/edit stack_interp.py
#Can run stack_interp on single node or parallel qsub ~/src/gmbtools/gmbtools/stack_interp_parallel.pbs
~/src/gmbtools/gmbtools/dem_post_mb_prep.sh

#Compute mass balance
#Update input filenames, output directory
~/src/gmbtools/gmbtools/mb_parallel.py

#Generate plots and run mb analysis
#~/src/gmbtools/gmbtools/mb_plot_gpd.py
Updated analysis in mb_plot_gpd.ipynb
mb_discharge.py - WBM analysis
mb_plot_gpd_compare.py - compare time periods, Brun and SRTM

#Mosaicking - run for both ASTER and WV/GE
#Need to update these steps

#1/3 arcsec vs 8 m AEA
#Add tile names
make_mos.sh
~/src/gmbtools/gmbtools/mos_eval.ipynb

#Create combined mosaic
make_mos_combined.sh
#Also review make_mos_filled.sh and site_rerun_demprep.sh for filling strategy 
