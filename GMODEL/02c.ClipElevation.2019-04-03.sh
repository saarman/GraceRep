#!/bin/bash
#SBATCH -p day 
##SBATCH -n 1 -c 1  -N 1  
##SBATCH -t 24:00:00
#SBATCH --output=02c.ClipElevation.2019-04-03.sh.log
#SBATCH --job-name=02c.ClipElevation.2019-04-03.sh
##DATA SOURCE: https://github.com/evlynpless/MOSQLAND/tree/master/MERIT/altitude
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/02c.ClipElevation.2019-04-03.sh
##LOG FILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/02c.ClipElevation.2019-04-03.sh.log
INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/MERIT/altitude
RASTERNAME=altitude_1KMmedian_MERIT

CLIPNAME=UgandaKenyaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin 28.6 4.73 42.5 -4.8 $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -999
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc

CLIPNAME=UgandaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin 28.6 4.73 35.4 -1.5 $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -999
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc

CLIPNAME=KenyaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin 33.7 4.73 42.5 -4.8 $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -999
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc
