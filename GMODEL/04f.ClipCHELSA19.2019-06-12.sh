#!/bin/bash
#SBATCH -p day 
##SBATCH -n 1 -c 1  -N 1  
##SBATCH -t 24:00:00
#SBATCH --output=04f.ClipCHELSA19.2019-06-12.sh.log
#SBATCH --job-name=04f.ClipCHELSA19.2019-06-12.sh
##DATA SOURCE: https://github.com/evlynpless/MOSQLAND/tree/master/bioclim
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/04f.ClipCHELSA19.2019-06-12.sh
##LOG FILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/04f.ClipCHELSA19.2019-06-12.sh.log
INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim
RASTERNAME=CHELSA_bio10

for n in $(seq 1 19 ); do
CLIPNAME=SouthAmericaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin -100.6 22.1 -34.47 -56.7 $INDIR/${RASTERNAME}_${n}.tif  $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif --msknodata -32768 -nodata -999   
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.asc
done
