#! /bin/bash

#Prepare to generate stacks for each RGI polygon using rgi_dem_trend.py
#Run after dem_post_parallel.pbs to align, then dem_align_post.py

#TODO
#Switch from shp to gpkg

#topdir=/nobackup/deshean/hma/dem_coreg
#topdir=/nobackup/deshean/hma/aster/dsm/dem_align_ASTERonly
topdir=/nobackup/deshean/hma/combined_aster_wv

if [ ! -d $topdir ] ; then 
    mkdir -pv $topdir
fi

cd $topdir

#prefix=dem_align_ASTER_round2
#prefix=dem_align_WV
prefix=dem_align_ASTER_WV
#prefix=dem_align_ASTERonly

#Throw out outliers
#mkdir dem_align_bad; for i in $(cat ${prefix}_bad_fn.txt); do mv -v $(echo $i | awk -F'/' '{print $1 "/" $2}') dem_align_bad/; done

#Organize into annual subdir, important for many files and shell list length limitations
valid_years=$(seq 2000 2018)

#ASTER
#valid_years=$(cat ${prefix}_good_fn.txt | cut -c 1-4 | sort -u)
#WV
#valid_years=$(cat ${prefix}_good_fn.txt | awk -F'/' '{print $2}' | cut -c 1-4 | sort -u)

yr1=$(echo $valid_years | awk '{print $1}')
yr2=$(echo $valid_years | awk '{print $NF}')
#yr_mid=2009

#ASTER round2
#for y in $valid_years ; do if [ ! -d years/$y ] ; then mkdir -pv years/$y ; fi ; cd years/$y ; for i in ../../${y}/*dem_align/*dem_align/*align.tif; do ln -s $i . ; done; cd ../../ ; done
#WV
#for y in $valid_years ; do if [ ! -d years/$y ] ; then mkdir -pv years/$y ; fi ; cd years/$y ; for i in ../../*track/${y}*dem_align/*align.tif; do ln -s $i . ; done; cd ../../ ; done

#Combined
for y in $valid_years
do 
    if [ ! -d years/$y ] ; then 
        mkdir -pv years/$y  
    fi 
    cd years/$y
    #WV
    for i in ../../../dem_coreg/*track/${y}*dem_align/*align.tif ; do ln -s $i . ; done
    #ASTER
    for i in ../../../aster/dsm/${y}/A*dem_align/*align.tif ; do ln -s $i . ; done
    #ASTER-only
    #for i in ../../../${y}/A*dem_align/*align.tif ; do ln -s $i . ; done
    cd ../../ 
done

parallel --progress "gdaltindex -t_srs EPSG:4326 years/{}/${prefix}_index_{}.shp years/{}/*align.tif" ::: $valid_years

if false ; then
    yr1=2007
    yr2=2018
    ogr_merge.sh ${prefix}_index_$yr1-$yr2.shp years/2*/${prefix}_index_*.shp
    shp_list=${prefix}_index_2007-2018.shp
else
    #Generate DEM index
    yr1=2000
    yr2=2018
    ogr_merge.sh ${prefix}_index_$yr1-$yr2.shp years/2*/${prefix}_index_*.shp

    yr1=2000
    yr2=2009
    ogr_merge.sh ${prefix}_index_$yr1-$yr2.shp years/2*/${prefix}_index_200[0-9].shp

    yr1=2009
    yr2=2018
    ogr_merge.sh ${prefix}_index_$yr1-$yr2.shp years/2*/${prefix}_index_2009.shp years/2*/${prefix}_index_201[0-9].shp

    #Now convert to aea projection
    shp_list="${prefix}_index_2000-2018.shp ${prefix}_index_2000-2009.shp ${prefix}_index_2009-2018.shp"
fi

parallel "ogr2ogr -t_srs '+proj=aea +lat_1=25 +lat_2=47 +lat_0=36 +lon_0=85 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ' {.}_aea.shp {}" ::: $shp_list 
shp_list=$(echo $shp_list | sed 's/.shp/_aea.shp/g')

#Prepare cmd files up front, in parallel
a1=0.0
a2=2.0
a3=9999.0
parallel --delay 10.0 "rgi_dem_trend.py {1} {2} {3}" ::: $shp_list ::: $a1 $a2 :::+ $a2 $a3

#Now run dem_post_parallel.pbs
