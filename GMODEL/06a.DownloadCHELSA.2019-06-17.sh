#!/bin/bash
#SBATCH -p day 
##SBATCH -n 1 -c 1  -N 1  
##SBATCH -t 24:00:00
#SBATCH --output=06a.DownloadCHELSA.2019-06-17.sh.log
#SBATCH --job-name=06a.DownloadCHELSA.2019-06-17.sh
##COMMAND: sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/06a.DownloadCHELSA.2019-06-17.sh
##LOG FILE: cat /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/06a.DownloadCHELSA.2019-06-17.sh.log
cd /home/fas/caccone/nps25/project/RASTERS/CHELSA_08-13

for VAR in `echo prec tmax tmean tmin`; do
	for y in $(seq 2008 2013 ); do
		for x in $(seq 1 9 ); do
			echo "https://www.wsl.ch/lud/chelsa/data/timeseries/${VAR}/CHELSA_${VAR}_${y}_0${x}_V1.2.1.tif"
			#wget https://www.wsl.ch/lud/chelsa/data/timeseries/${VAR}/CHELSA_${VAR}_${y}_0${x}_V1.2.1.tif
		done
		for n in $(seq 10 12 ); do
			echo "https://www.wsl.ch/lud/chelsa/data/timeseries/${VAR}/CHELSA_${VAR}_${y}_${x}_V1.2.1.tif"
			#wget https://www.wsl.ch/lud/chelsa/data/timeseries/${VAR}/CHELSA_${VAR}_${y}_${x}_V1.2.1.tif
		done
	done
done

