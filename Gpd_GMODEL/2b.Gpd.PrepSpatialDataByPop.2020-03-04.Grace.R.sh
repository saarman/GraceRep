#!/bin/bash
#SBATCH -p day 
#SBATCH -n 1 -c 4
#SBATCH --mem-per-cpu=3G
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/loomis/scratch60/fas/caccone/nps25/stdout/2b.Gpd.PrepSpatialDataByPop.2020-03-04.Grace.R.sh.%J.out
#SBATCH -o /gpfs/loomis/scratch60/fas/caccone/nps25/stderr/2b.Gpd.PrepSpatialDataByPop.2020-03-04.Grace.R.sh.%J.err
#SBATCH --output=2b.Gpd.PrepSpatialDataByPop.2020-03-04.Grace.R.sh.log
#SBATCH --job-name=2b.Gpd.PrepSpatialDataByPop.2020-03-04.Grace.R.sh
#sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/Gpd_GMODEL/2b.Gpd.PrepSpatialDataByPop.2020-03-04.Grace.R.sh

#load conda environment:
module purge
module load miniconda
source activate parallel_r

#run R script:
cd /home/fas/caccone/nps25/scripts/GraceRep/Gpd_GMODEL
R --vanilla --no-readline -q -f 2b.Gpd.PrepSpatialDataByPop.2020-03-04.Grace.R
###############################################
