
#Based on Evie's code at https://github.com/evlynpless/MOSQLAND/blob/master/RF/sc06_iterativeRF_Florida.R

library("sp")
library("spatstat")
library("maptools")
library("raster")
library("randomForest")
library("gdistance")
library("SDraw")
library("tidyverse")

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system

###############################################
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

G.table <- read.table(file="/project/fas/caccone/nps25/GMODEL/dist_matrix.csv", sep=",", header=T) # ... load coordinates and Fst, etc from file

begin.table <- G.table[,c(5,3)]
begin.coord <- begin.table
coordinates(begin.coord) <- c("long1", "lat1")

end.table <- G.table[,c(6,4)]
end.coord <- end.table
coordinates(end.coord) <- c("long2", "lat2")

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

spatial.p <- as(p, "SpatialLines")
proj4string(spatial.p) <- crs.geo  # define projection system of our data

print("spatial lines done")

###############################################
#Create raster stack 
###############################################

raster01I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/pet_mean_UgandaClip.tif")
raster02I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/slope_1KMmedian_MERIT_UgandaClip.tif")
raster03I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/AI_annual_UgandaClip.tif")
raster04I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/altitude_1KMmedian_MERIT_UgandaClip.tif")
raster05I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_10_UgandaClip.tif")
raster06I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_11_UgandaClip.tif")
raster07I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_12_UgandaClip.tif")
raster08I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_13_UgandaClip.tif")
raster09I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_14_UgandaClip.tif")
raster10I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_15_UgandaClip.tif")
raster11I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_16_UgandaClip.tif")
raster12I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_17_UgandaClip.tif")
raster13I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_18_UgandaClip.tif")
raster14I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_19_UgandaClip.tif")
raster15I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_1_UgandaClip.tif")
raster16I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_2_UgandaClip.tif")
raster17I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_3_UgandaClip.tif")
raster18I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_4_UgandaClip.tif")
raster19I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_5_UgandaClip.tif")
raster20I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_6_UgandaClip.tif")
raster21I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_7_UgandaClip.tif")
raster22I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_8_UgandaClip.tif")
raster23I = raster("/project/fas/caccone/nps25/RASTERS/UgandaClips/CHELSA_bio10_9_UgandaClip.tif")

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
          raster23)

names(env) [1]  <- "pet_mean"
names(env) [2]  <- "slope"
names(env) [3]  <- "AI_annual"
names(env) [4]  <- "altitude"
names(env) [5]  <- "CHELSA_bio10_10"
names(env) [6]  <- "CHELSA_bio10_11"
names(env) [7]  <- "CHELSA_bio10_12"
names(env) [8]  <- "CHELSA_bio10_13"
names(env) [9]  <- "CHELSA_bio10_14"
names(env) [10] <- "CHELSA_bio10_15"
names(env) [11] <- "CHELSA_bio10_16"
names(env) [12] <- "CHELSA_bio10_17"
names(env) [13] <- "CHELSA_bio10_18"
names(env) [14] <- "CHELSA_bio10_19"
names(env) [15] <- "CHELSA_bio10_1"
names(env) [16] <- "CHELSA_bio10_2"
names(env) [17] <- "CHELSA_bio10_3"
names(env) [18] <- "CHELSA_bio10_4"
names(env) [19] <- "CHELSA_bio10_5"
names(env) [20] <- "CHELSA_bio10_6"
names(env) [21] <- "CHELSA_bio10_7"
names(env) [22] <- "CHELSA_bio10_8"
names(env) [23] <- "CHELSA_bio10_9"

print("raster stack done")

########################################
#Calculate mean of straight lines and making initial RF model
#######################################
StraightMean <- raster::extract(env, spatial.p, fun=mean, na.rm=TRUE)

StraightMeanDF <- as.data.frame(StraightMean)

StraightMeanDF$FST_arl <- G.table$FST_arl

#option of trying DPS
#StraightMeanDF$DPS <- G.table$DPS

Straight_RF1 = randomForest(FST_arl ~ pet_mean +
                                      slope +
                                      AI_annual +
                                      altitude +
                                      CHELSA_bio10_10 +
                                      CHELSA_bio10_11 +
                                      CHELSA_bio10_12 +
                                      CHELSA_bio10_13 +
                                      CHELSA_bio10_14 +
                                      CHELSA_bio10_15 +
                                      CHELSA_bio10_16 +
                                      CHELSA_bio10_17 +
                                      CHELSA_bio10_18 +
                                      CHELSA_bio10_19 +
                                      CHELSA_bio10_1 +
                                      CHELSA_bio10_2 +
                                      CHELSA_bio10_3 +
                                      CHELSA_bio10_4 +
                                      CHELSA_bio10_5 +
                                      CHELSA_bio10_6 +
                                      CHELSA_bio10_7 +
                                      CHELSA_bio10_8 +
                                      CHELSA_bio10_9, importance=TRUE, na.action=na.omit, data=StraightMeanDF)

Straight_RF

StraightPred <- predict(env, Straight_RF)

print("first prediction resistance surface done")

pred.cond <- 1/StraightPred #build conductance surface

save.image(file = "/project/fas/caccone/nps25/GMODEL/03a.GeneticRandomForest.2019-04-06.RData")

test = summary(LCP_RF)
