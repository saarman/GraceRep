#!/bin/bash
#SBATCH -p day 
##SBATCH -n 1 -c 1  -N 1  
##SBATCH -t 24:00:00
#SBATCH --output=04e.ClipPET.2019-06-12.sh.log
#SBATCH --job-name=04e.ClipPET.2019-06-12.sh
##DATA SOURCE: https://github.com/evlynpless/MOSQLAND/tree/master/PET
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/04e.ClipPET.2019-06-12.sh
##LOG FILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/04e.ClipPET.2019-06-12.sh.log
INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/PET
RASTERNAME=pet_mean

CLIPNAME=SouthAmericaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin -100.6 22.1 -34.47 -56.7 $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif --msknodata -16257 -nodata -999
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc
