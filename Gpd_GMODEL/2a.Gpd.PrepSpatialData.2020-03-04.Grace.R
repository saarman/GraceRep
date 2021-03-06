###############################################
#R Script: Note that this is also tested locally as "~/Dropbox/Caccone_Aksoy/Glossina-spatial/Gpd-RF-Genetic-Model/2a.Gpd.PrepSpatialData.2020-03-04.Rmd"
###############################################
#first thing you do is load the library that is for this analysis
library("raster")
library("rgdal")
library("dismo")
library("XML")
library("maps")
library("sp")
library("spatstat")
library("maptools")
library("randomForest")
library("gdistance")
library("dplyr")
library("foreach")
library("doParallel")
library("doMC")

#set working directory
setwd("~/project/Gpd_GMODEL")

#set coordinate system for mapping
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 

###############################################>
#Create raster stack 
###############################################>
#Load a stack of raster layers
rast_stack <- stack("/home/nps25/project/RASTERS/KenyaClips/chelsa_merit_vars_kenya.tif") 
names(rast_stack) <- c(paste0("BIO",c(8:11),"S"),paste0("BIO",c(16:19),"S"),paste0("BIO",c(1:19)),"slope","altitude")

#isolate each raster layer
###############################################>
raster01I = subset(rast_stack, 1)
raster02I = subset(rast_stack, 2)
raster03I = subset(rast_stack, 3)
raster04I = subset(rast_stack, 4)
raster05I = subset(rast_stack, 5)
raster06I = subset(rast_stack, 6)
raster07I = subset(rast_stack, 7)
raster08I = subset(rast_stack, 8)
raster09I = subset(rast_stack, 9)
raster10I = subset(rast_stack, 10)
raster11I = subset(rast_stack, 11)
raster12I = subset(rast_stack, 12)
raster13I = subset(rast_stack, 13)
raster14I = subset(rast_stack, 14)
raster15I = subset(rast_stack, 15)
raster16I = subset(rast_stack, 16)
raster17I = subset(rast_stack, 17)
raster18I = subset(rast_stack, 18)
raster19I = subset(rast_stack, 19)
raster20I = subset(rast_stack, 20)
raster21I = subset(rast_stack, 21)
raster22I = subset(rast_stack, 22)
raster23I = subset(rast_stack, 23)
raster24I = subset(rast_stack, 24)
raster25I = subset(rast_stack, 25)
raster26I = subset(rast_stack, 26)
raster27I = subset(rast_stack, 27)
raster28I = subset(rast_stack, 28)
raster29I = subset(rast_stack, 29)
raster30I = (raster29I*0)+1 #all pixels are 1

#The "I" times 1 puts the raster in active memory
###############################################>
raster01 <- raster01I*1
raster02 <- raster02I*1
raster03 <- raster03I*1
raster04 <- raster04I*1
raster05 <- raster05I*1
raster06 <- raster06I*1
raster07 <- raster07I*1
raster08 <- raster08I*1
raster09 <- raster09I*1
raster10 <- raster10I*1
raster11 <- raster11I*1
raster12 <- raster12I*1
raster13 <- raster13I*1
raster14 <- raster14I*1
raster15 <- raster15I*1
raster16 <- raster16I*1
raster17 <- raster17I*1
raster18 <- raster18I*1
raster19 <- raster19I*1
raster20 <- raster20I*1
raster21 <- raster21I*1
raster22 <- raster22I*1
raster23 <- raster23I*1
raster24 <- raster24I*1
raster25 <- raster25I*1
raster26 <- raster26I*1
raster27 <- raster27I*1
raster28 <- raster28I*1
raster29 <- raster29I*1
raster30 <- raster30I*1

#make sure the projection is correct for each:
###############################################>
proj4string(raster01) <- crs.geo
proj4string(raster02) <- crs.geo
proj4string(raster03) <- crs.geo
proj4string(raster04) <- crs.geo
proj4string(raster05) <- crs.geo
proj4string(raster06) <- crs.geo
proj4string(raster07) <- crs.geo
proj4string(raster08) <- crs.geo
proj4string(raster09) <- crs.geo
proj4string(raster10) <- crs.geo
proj4string(raster11) <- crs.geo
proj4string(raster12) <- crs.geo
proj4string(raster13) <- crs.geo
proj4string(raster14) <- crs.geo
proj4string(raster15) <- crs.geo
proj4string(raster16) <- crs.geo
proj4string(raster17) <- crs.geo
proj4string(raster18) <- crs.geo
proj4string(raster19) <- crs.geo
proj4string(raster20) <- crs.geo
proj4string(raster21) <- crs.geo
proj4string(raster22) <- crs.geo
proj4string(raster23) <- crs.geo
proj4string(raster24) <- crs.geo
proj4string(raster25) <- crs.geo
proj4string(raster26) <- crs.geo
proj4string(raster27) <- crs.geo
proj4string(raster28) <- crs.geo
proj4string(raster29) <- crs.geo
proj4string(raster30) <- crs.geo

###############################################>
env=stack(raster01,
          raster02,
          raster03,
          raster04,
          raster05,
          raster06,
          raster07,
          raster08,
          raster09,
          raster10,
          raster11,
          raster12,
          raster13,
          raster14,
          raster15,
          raster16,
          raster17,
          raster18,
          raster19,
          raster20,
          raster21,
          raster22,
          raster23,
          raster24,
          raster25,
          raster26,
          raster27,
          raster28,
          raster29)
pixels <- raster30
###############################################>

##BIO1-BIO19 are the bioclim variables. If it has an "S" at the end, it is the seasonal (rather than quarterly) adaptation of the corresponding bioclim variable.
names(env) [1] <-   "Bio08S_temp_mean_wettestS"
names(env) [2] <-   "Bio09S_temp_mean_driestS"  
names(env) [3] <-   "Bio10S_temp_mean_warmestS"
names(env) [4] <-   "Bio11S_temp_mean_coldestS"
names(env) [5] <-   "Bio16S_precip_wettestS"
names(env) [6] <-   "Bio17S_precip_driestS"
names(env) [7] <-   "Bio18S_precip_warmestS"
names(env) [8] <-   "Bio19S_precip_coldestS"
names(env) [9]  <-  "Bio01_temp_mean_annual"
names(env) [10]  <- "Bio02_temp_mean_diurnal_range"
names(env) [11]  <- "Bio03_isothermality"
names(env) [12]  <- "Bio04_temp_seasonality"
names(env) [13]  <- "Bio05_temp_max_warmestMo"
names(env) [14]  <- "Bio06_temp_min_coldestMo"
names(env) [15]  <- "Bio07_temp_annual_range"
names(env) [16]  <- "Bio08_temp_mean_wettestQ"
names(env) [17]  <- "Bio09_temp_mean_driestQ"  
names(env) [18]  <- "Bio10_temp_mean_warmestQ"
names(env) [19]  <- "Bio11_temp_mean_coldestQ"
names(env) [20]  <- "Bio12_precip_annual"
names(env) [21]  <- "Bio13_precip_wettestMo"
names(env) [22]  <- "Bio14_precip_driestMo"
names(env) [23]  <- "Bio15_precip_seasonality"
names(env) [24]  <- "Bio16_precip_wettestQ"
names(env) [25]  <- "Bio17_precip_driestQ"
names(env) [26]  <- "Bio18_precip_warmestQ"
names(env) [27]  <- "Bio19_precip_coldestQ"
names(env) [28]  <- "slope"
names(env) [29]  <- "altitude"
names(pixels) <- "pixels"

#NOTE: If you add future, the names of the rasters in the stack have to be the same as the current ones used in building the model.
#NOTE: current temperature rasters are in ˚K*10. future temperature rasters are in ˚C*10... so the conversion from future to current for temperature is current=raster01F/10)+273.15)*10.
###############################################>

# clean up objects
x <- c("raster01I")
x <- c(x, "raster02I")
x <- c(x, "raster03I")
x <- c(x, "raster04I")
x <- c(x, "raster05I")
x <- c(x, "raster06I")
x <- c(x, "raster07I")
x <- c(x, "raster08I")
x <- c(x, "raster09I")
x <- c(x, "raster10I")
x <- c(x, "raster11I")
x <- c(x, "raster12I")
x <- c(x, "raster13I")
x <- c(x, "raster14I")
x <- c(x, "raster15I")
x <- c(x, "raster16I")
x <- c(x, "raster17I")
x <- c(x, "raster18I")
x <- c(x, "raster19I")
x <- c(x, "raster20I")
x <- c(x, "raster21I")
x <- c(x, "raster22I")
x <- c(x, "raster23I")
x <- c(x, "raster24I")
x <- c(x, "raster25I")
x <- c(x, "raster26I")
x <- c(x, "raster27I")
x <- c(x, "raster28I")
x <- c(x, "raster29I")
x <- c(x, "raster30I")
x <- c(x, "raster01")
x <- c(x, "raster02")
x <- c(x, "raster03")
x <- c(x, "raster04")
x <- c(x, "raster05")
x <- c(x, "raster06")
x <- c(x, "raster07")
x <- c(x, "raster08")
x <- c(x, "raster09")
x <- c(x, "raster10")
x <- c(x, "raster11")
x <- c(x, "raster12")
x <- c(x, "raster13")
x <- c(x, "raster14")
x <- c(x, "raster15")
x <- c(x, "raster16")
x <- c(x, "raster17")
x <- c(x, "raster18")
x <- c(x, "raster19")
x <- c(x, "raster20")
x <- c(x, "raster21")
x <- c(x, "raster22")
x <- c(x, "raster23")
x <- c(x, "raster24")
x <- c(x, "raster25")
x <- c(x, "raster26")
x <- c(x, "raster27")
x <- c(x, "raster28")
x <- c(x, "raster29")
x <- c(x, "raster30")
rm(list=x)

# map a few variables
###############################################>
plot(env$Bio19S_precip_coldestS,axes=FALSE, box=FALSE)
map(,ylim=c(-4.8, 5),xlim=c(33.7, 42.5),add=T)
plot(pixels,axes=FALSE, box=FALSE)
map(,ylim=c(-4.8, 5),xlim=c(33.7, 42.5),add=T)

###############################################>
# individual pairwise data
###############################################>
# load individual pairwise data and calculate spatial lines
G.table.all <- read.table(file="Gpd_KenTza_11loci_659indv_indiv_pairwise.txt",header=T)
G.table <- data.frame(G.table.all[G.table.all$Ind1_siteID!=G.table.all$Ind2_siteID & G.table.all$Ind1_cluster==G.table.all$Ind2_cluster,])
unique.G.table <- unique.data.frame(G.table.all[G.table.all$Ind1_siteID!=G.table.all$Ind2_siteID & G.table.all$Ind1_cluster==G.table.all$Ind2_cluster,c("Ind1_pixelLong","Ind1_pixelLat","Ind2_pixelLong","Ind2_pixelLat")])

###############################################>
begin.table <- unique.G.table[,c("Ind1_pixelLong","Ind1_pixelLat")]
begin.coord <- begin.table
begin.Sp.points <- SpatialPoints(begin.coord)

end.table <- unique.G.table[,c("Ind2_pixelLong","Ind2_pixelLat")]
end.coord <- end.table
end.Sp.points <- SpatialPoints(end.coord)
###############################################>
l <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))
spatial.l <- as(l, "SpatialLines")
proj4string(spatial.l) <- crs.geo  # define projection system of our data

# Plot unique lines as SpatialLines
###############################################>
plot(env$Bio19S_precip_coldestS,axes=FALSE, box=FALSE)
map(,ylim=c(-4.8, 5),xlim=c(33.7, 42.5),add=T)
plot(spatial.l,add=T)

# Calculate mean of unique lines for each variable
registerDoMC(cores=detectCores()) #get parallelization set up
StraightMean <- raster::extract(env, spatial.l, fun=mean, na.rm=TRUE)
summary(StraightMean)
gc()

# calculate sum of pixels along straight lines
###############################################>
StraightSumPixels <- raster::extract(pixels, spatial.l, fun=sum, na.rm=TRUE)
summary(StraightSumPixels) 

# create a final data table
###############################################>
StriaghtLinesTable <- cbind(unique.G.table,as.data.frame(StraightMean),as.data.frame(StraightSumPixels))
StraightLinesDF <- left_join(G.table,StriaghtLinesTable)

save.image(file="2a.Gpd.PrepSpatialData.2020-03-04.Grace.RData")
