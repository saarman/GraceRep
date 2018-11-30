#!/bin/bash
#SBATCH -p day 
#SBATCH -n 1 -c 1  -N 1  
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/loomis/scratch60/fas/caccone/nps25/stdout/01a.CreateSpatialLines.2018-11-30.R.sh.%J.out
#SBATCH -o /gpfs/loomis/scratch60/fas/caccone/nps25/stderr/01a.CreateSpatialLines.2018-11-30.R.sh.%J.err
#SBATCH --output=01a.CreateSpatialLines.2018-11-30.R.sh.log
#SBATCH --job-name=01a.CreateSpatialLines.2018-11-30.R.sh
#sbatch /gpfs/loomis/home.grace/fas/caccone/nps25/scripts/GraceRep/GMODEL/01a.CreateSpatialLines.2018-11-30.R.sh

#export MERIT=/project/fas/sbsc/ga254/grace0.grace.hpc.yale.internal/dataproces/MERIT
#export SCRATCH=/gpfs/scratch60/fas/sbsc/ga254/grace0/dataproces/MERIT_BK
#export RAM=/dev/shm
#export KM=5.00

module load Apps/R/3.3.2-generic

R --vanilla --no-readline   -q  <<'EOF'

#R Script: Note that this is also tested and saved locally as "/Users/Norah/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-PopGen/LoadSpatialLines.R"
###############################################
#source(01a.CreateSpatialLines.2018-11-30.R)

#First install packages in interactive mode, then load packages:
library("sp")
library("spatstat")
library("maptools")

#set working directory
setwd("/home/fas/caccone/nps25/project/GMODEL")

#One way to do it is to use the coordinates function:
G.table <- read.table(file="dist_matrix.csv", sep=",", header=T) #-> ... load coordinates and Fst, etc from file 
coordinates(G.table) <- c("long1", "lat1") #->
coordinates1 <- coordinates(G.table) #->
plot(coordinates1)

#Another way to do it is with SpatialPoints
G.table <- read.table(file="dist_matrix.csv", sep=",", header=T)  #->

G.coordinates1 <- G.table[,c(5,3)] #-> ... select the correct collumns 
G.points1 <- SpatialPoints(G.coordinates1)  #-> ... converted into a spatial object
plot(G.points1)

summary(G.points1)
is.projected(G.points1)
crs.geo <- CRS("+init=EPSG:32633 +datum=WGS84") #-> ... add coordinate system
crs.geo
proj4string(G.points1) <- crs.geo  #-> define projection system of our data
summary(G.points1)
is.projected(G.points1)

G.coordinates2 <- G.table[,c(6,4)] #->
G.points2 <- SpatialPoints(G.coordinates2) #->
plot(G.points2)

summary(G.points2)
is.projected(G.points2)
crs.geo <- CRS("+init=EPSG:32633 +datum=WGS84") #-> ... add coordinate system
proj4string(G.points2) <- crs.geo  #-> define projection system of our data
summary(G.points2)
is.projected(G.points2)

#Now a way to create lines from table --
#create dataframes of begin and end coordinates:
G.table <- read.table(file="dist_matrix.csv", sep=",", header=T) #-> ... load coordinates and Fst, etc from file

begin.table <- G.table[,c(5,3)] #->
begin.coord <- begin.table #->
coordinates(begin.coord) <- c("long1", "lat1") #->

end.table <- G.table[,c(6,4)] #->
end.coord <- end.table #->
coordinates(end.coord) <- c("long2", "lat2") #->

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2])))) #->
plot(p)

spatial.p <- as(p, "SpatialLines") #->
crs.geo <- CRS("+init=EPSG:32633 +datum=WGS84") #-> ... add coordinate system
proj4string(spatial.p) <- crs.geo  #-> define projection system of our data
summary(spatial.p)

plot(spatial.p)

save.image("01a.CreateSpatialLines.2018-11-30.Rdata")

#To load the saved image later:
# library("sp")
# library("spatstat")
# library("maptools")
# setwd("/home/fas/caccone/nps25/project/GMODEL")
# load("01a.CreateSpatialLines.2018-11-30.Rdata")
# summary(spatial.p)
# plot(p)
# plot(spatial.p)
###############################################
EOF
