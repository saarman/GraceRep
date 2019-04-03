#!/bin/bash
#SBATCH -p day 
#SBATCH -n 1 -c 1  -N 1  
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/loomis/scratch60/fas/caccone/nps25/stdout/02a.ClipAridity.2019-04-03.sh.out
#SBATCH -o /gpfs/loomis/scratch60/fas/caccone/nps25/stderr/02a.ClipAridity.2019-04-03.sh.err
#SBATCH --output=02a.ClipAridity.2019-04-03.sh.log
#SBATCH --job-name=02a.ClipAridity.2019-04-03.sh
## sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/02a.ClipAridity.2019-04-03.sh
INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY
RASTERNAME=AI_annual

CLIPNAME=UgandaKenyaClip
OUTDIR=/project/caccone/nps25/project/RASTERS/${CLIPNAME}s
gdal_translate  -projwin -4.8 42.5 4.73 28.6  $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i ${RASTERNAME}_${CLIPNAME}.tif -m ${RASTERNAME}_${CLIPNAME}.tif -o ${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -9999   
gdal_translate -of AAIGrid $OUTDIR/A${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc

CLIPNAME=UgandaClip
OUTDIR=/project/caccone/nps25/project/RASTERS/${CLIPNAME}s
gdal_translate  -projwin -1.5 35.4 4.73 28.6  $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i ${RASTERNAME}_${CLIPNAME}.tif -m ${RASTERNAME}_${CLIPNAME}.tif -o ${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -9999   
gdal_translate -of AAIGrid $OUTDIR/A${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc

CLIPNAME=KenyaClip
OUTDIR=/project/caccone/nps25/project/RASTERS/${CLIPNAME}s
gdal_translate  -projwin -4.8 42.5 4.73 33.7  $INDIR/${RASTERNAME}.tif  $OUTDIR/${RASTERNAME}_${CLIPNAME}.tif
pksetmask -i ${RASTERNAME}_${CLIPNAME}.tif -m ${RASTERNAME}_${CLIPNAME}.tif -o ${RASTERNAME}_${CLIPNAME}.tif --msknodata -1 -nodata -9999   
gdal_translate -of AAIGrid $OUTDIR/A${RASTERNAME}_${CLIPNAME}.tif $OUTDIR/${RASTERNAME}_${CLIPNAME}.asc
