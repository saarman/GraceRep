#!/bin/bash
#SBATCH -p day 
##SBATCH -n 1 -c 1  -N 1  
##SBATCH -t 24:00:00
#SBATCH --output=04a.ClipAridity.2019-12-10.sh.log
#SBATCH --job-name=04a.ClipAridity.2019-12-10.sh
##DATA SOURCE: https://github.com/evlynpless/MOSQLAND/tree/master/ARIDITY
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/04a.ClipAridity.2019-12-10.sh
##LOG FILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/04a.ClipAridity.2019-12-10.sh.log
module load PKTOOLS/2.6.7.6-foss-2018a
module load GDAL/2.2.3-foss-2018a-Python-2.7.14
INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY
RASTERNAME=AI_annual
CLIPNAME=CaliforniaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin -125 43 -113 32 $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -999
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc
