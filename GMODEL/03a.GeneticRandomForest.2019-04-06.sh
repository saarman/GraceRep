#!/bin/bash
#SBATCH -p day 
#SBATCH -n 1 -c 1  -N 1  
#SBATCH -t 24:00:00
#SBATCH --output=03a.GeneticRandomForest.2019-04-06.sh.log
#SBATCH --job-name=03a.GeneticRandomForest.2019-04-06.sh
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/03a.GeneticRandomForest.2019-04-06.sh
##LOGFILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/03a.GeneticRandomForest.2019-04-06.sh.log

module load Apps/R/3.3.2-generic
module load Rpkgs/RGDAL/1.2-5
R --vanilla -no-readline -q  -f  /home/fas/caccone/nps25/scripts/GraceRep/GMODEL/03a.GeneticRandomForest.2019-04-06.R
