
#R Script: Note that this is also tested and saved locally as "/Users/Norah/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-PopGen/LoadSpatialLines.R"
###############################################
#First install packages: These only seem to install in temp file, can we get them to last?
install.packages("sp", repos="https://mirrors.sorengard.com/cran/")
install.packages("spatstat", repos="https://mirrors.sorengard.com/cran/")

#now load packages
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

#Need to save the output in some way since this will not be interactive in final form
