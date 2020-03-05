#!/bin/bash
#SBATCH -p bigmem
#SBATCH --mem=500g
#SBATCH -n 1 -c 8 -N 1
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anusha.bishop@yale.edu
#SBATCH --job-name=sc02_pt2_dpsKEN_LCP.sh
#SBATCH --output=sc02_pt2_dpsKEN_LCP.sh.log

####  for foldnum in 1 2 3 4 5 6 7 8 9 10; do sbatch --export=foldnum=$foldnum /home/fas/caccone/nps25/scripts/GraceRep/Gpd_GMODEL/sc02_pt2_dpsKEN_LCP_nps.sh  ; done

#for testing script
####  for foldnum in 1; do sbatch --export=foldnum=$foldnum /home/fas/caccone/apb56/scripts/GPDGENCON/sc02_pt2_dpsgeoKEN_LCP.sh  ; done

ulimit -c 0

module purge
module load miniconda
source activate parallel_r

# --slave      use if you only want to see output

export foldnum=$foldnum

#run script with just envvars
R --vanilla --no-readline -q  -f /home/fas/caccone/nps25/scripts/GraceRep/Gpd_GMODEL/sc02_pt2_dpsgeoKEN_LCP_nps.R
