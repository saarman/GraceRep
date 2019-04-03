#!/bin/bash
#SBATCH -p day 
##SBATCH -n 1 -c 1  -N 1  
##SBATCH -t 24:00:00
#SBATCH --output=02f.ClipCHELSA19.2019-04-03.sh.log
#SBATCH --job-name=02f.ClipCHELSA19.2019-04-03.sh
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/02f.ClipCHELSA19.2019-04-03.sh
##LOG FILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/02f.ClipCHELSA19.2019-04-03.sh.log
INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim
RASTERNAME=CHELSA_bio10

for n in $(seq 1 19 ); do
CLIPNAME=UgandaKenyaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin 28.6 4.73 42.5 -4.8 $INDIR/${RASTERNAME}_${n}.tif  $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif --msknodata -32768 -nodata -999  
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.asc

CLIPNAME=UgandaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin 28.6 4.73 35.4 -1.5 $INDIR/${RASTERNAME}_${n}.tif  $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif --msknodata -32768 -nodata -999   
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.asc

CLIPNAME=KenyaClip
OUTDIR=/project/fas/caccone/nps25/RASTERS/${CLIPNAME}s
gdal_translate  -projwin 33.7 4.73 42.5 -4.8 $INDIR/${RASTERNAME}_${n}.tif  $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif
pksetmask -i $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -m $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif -o $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif --msknodata -32768 -nodata -999   
gdal_translate -of AAIGrid $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${n}_${CLIPNAME}.asc
done
