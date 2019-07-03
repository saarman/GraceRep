#Based on Evie's code at https://github.com/evlynpless/MOSQLAND/blob/master/RF/sc06_iterativeRF_Florida.R

library("classInt")
library("raster")
library("rgdal")
library("dismo")
library("XML")
library("maps")
library("sp")
library("classInt")

library("sp")
library("spatstat")
library("maptools")
library("raster")
library("randomForest")
library("gdistance")
library("SDraw")
library("tidyverse")

setwd("~/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-RF-Genetic-Model")

crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # ... add coordinate system

###############################################
#Plot lines as SpatialLines:
###############################################

#Plot straigt lines for first iteration of RF

G.table <- read.table(file="dist_matrix.csv", sep=",", header=T) # ... load coordinates and Fst, etc from file

begin.table <- G.table[,c(5,3)]
begin.coord <- begin.table
begin.Sp.points <- SpatialPoints(begin.coord)

end.table <- G.table[,c(6,4)]
end.coord <- end.table
end.Sp.points <- SpatialPoints(end.coord)

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))

plot(p)

spatial.p <- as(p, "SpatialLines")
proj4string(spatial.p) <- crs.geo  # define projection system of our data

print("spatial lines done")

###############################################
#Calculate geodist for each line and add it to the G.table
###############################################

geodist <- pointDistance(begin.Sp.points, end.Sp.points, lonlat=TRUE, allpairs=FALSE)
summary(geodist)
head(geodist)

G.table <- cbind(G.table,geodist)


###############################################
#Create raster stack 
###############################################
#The stack of bioclim rasters from Anusha is in order 1-19
#The stack of custom seasonal biovars are in this order: equivalent of 8-11 and 16-19 bioclim quarters

#Current Original Bioclim Variables: /home/fas/caccone/apb56/project/CHELSA/biovars/UgandaKenyaBiovarsAllYears.tif
UgandaKenyaBiovarsAllYears = stack("./Uganda08-13/UgandaKenyaBiovarsAllYears.tif")
UgandaBiovarsAllYears <- crop(UgandaKenyaBiovarsAllYears, extent(28.6, 35.4, -1.5, 4.73))

#Current Seasonal (New) Biovars: /home/fas/caccone/apb56/project/CHELSA/biovars/UgandaKenyaBiovarsSeasonalAllYears.tif
UgandaKenyaBiovarsSeasonalAllYears = stack("./Uganda08-13/UgandaKenyaBiovarsSeasonalAllYears.tif")
UgandaBiovarsSeasonalAllYears <- crop(UgandaKenyaBiovarsSeasonalAllYears, extent(28.6, 35.4, -1.5, 4.73))

#Future Original Bioclim Variables: /home/fas/caccone/apb56/project/CHELSA/biovars/UgandaBiovarsFuture.tif
UgandaOrigBiovarsFuture = stack("./UgandaRCP4.5/UgandaBiovarsFuture.tif")
UgandaBiovarsFuture <- crop(UgandaOrigBiovarsFuture, extent(28.6, 35.4, -1.5, 4.73))

#Future Seasonal (New) Biovars: /home/fas/caccone/apb56/project/CHELSA/biovars/UgandaBiovarsSeasonalFuture.tif
UgandaOrigBiovarsSeasonalFuture = stack("./UgandaRCP4.5/UgandaBiovarsSeasonalFuture.tif")
UgandaBiovarsSeasonalFuture <- crop(UgandaOrigBiovarsSeasonalFuture, extent(28.6, 35.4, -1.5, 4.73))

raster01I = subset(UgandaBiovarsAllYears, 1)
raster02I = subset(UgandaBiovarsAllYears, 2)
raster03I = subset(UgandaBiovarsAllYears, 3)
raster04I = subset(UgandaBiovarsAllYears, 4)
raster05I = subset(UgandaBiovarsAllYears, 5)
raster06I = subset(UgandaBiovarsAllYears, 6)
raster07I = subset(UgandaBiovarsAllYears, 7)
raster08I = subset(UgandaBiovarsAllYears, 8)
raster09I = subset(UgandaBiovarsAllYears, 9)
raster10I = subset(UgandaBiovarsAllYears, 10)
raster11I = subset(UgandaBiovarsAllYears, 11)
raster12I = subset(UgandaBiovarsAllYears, 12)
raster13I = subset(UgandaBiovarsAllYears, 13)
raster14I = subset(UgandaBiovarsAllYears, 14)
raster15I = subset(UgandaBiovarsAllYears, 15)
raster16I = subset(UgandaBiovarsAllYears, 16)
raster17I = subset(UgandaBiovarsAllYears, 17)
raster18I = subset(UgandaBiovarsAllYears, 18)
raster19I = subset(UgandaBiovarsAllYears, 19)
raster20I = subset(UgandaBiovarsSeasonalAllYears, 1)
raster21I = subset(UgandaBiovarsSeasonalAllYears, 2)
raster22I = subset(UgandaBiovarsSeasonalAllYears, 3)
raster23I = subset(UgandaBiovarsSeasonalAllYears, 4)
raster24I = subset(UgandaBiovarsSeasonalAllYears, 5)
raster25I = subset(UgandaBiovarsSeasonalAllYears, 6)
raster26I = subset(UgandaBiovarsSeasonalAllYears, 7)
raster27I = subset(UgandaBiovarsSeasonalAllYears, 8)
raster28I = raster("./Uganda08-13/altitude_1KMmedian_MERIT_UgandaClip.tif")
raster29I = raster("./Uganda08-13/slope_1KMmedian_MERIT_UgandaClip.tif")

#and future:
raster01FI = subset(UgandaBiovarsFuture, 1)
raster02FI = subset(UgandaBiovarsFuture, 2)
raster03FI = subset(UgandaBiovarsFuture, 3)
raster04FI = subset(UgandaBiovarsFuture, 4)
raster05FI = subset(UgandaBiovarsFuture, 5)
raster06FI = subset(UgandaBiovarsFuture, 6)
raster07FI = subset(UgandaBiovarsFuture, 7)
raster08FI = subset(UgandaBiovarsFuture, 8)
raster09FI = subset(UgandaBiovarsFuture, 9)
raster10FI = subset(UgandaBiovarsFuture, 10)
raster11FI = subset(UgandaBiovarsFuture, 11)
raster12FI = subset(UgandaBiovarsFuture, 12)
raster13FI = subset(UgandaBiovarsFuture, 13)
raster14FI = subset(UgandaBiovarsFuture, 14)
raster15FI = subset(UgandaBiovarsFuture, 15)
raster16FI = subset(UgandaBiovarsFuture, 16)
raster17FI = subset(UgandaBiovarsFuture, 17)
raster18FI = subset(UgandaBiovarsFuture, 18)
raster19FI = subset(UgandaBiovarsFuture, 19)
raster20FI = subset(UgandaBiovarsSeasonalFuture, 1)
raster21FI = subset(UgandaBiovarsSeasonalFuture, 2)
raster22FI = subset(UgandaBiovarsSeasonalFuture, 3)
raster23FI = subset(UgandaBiovarsSeasonalFuture, 4)
raster24FI = subset(UgandaBiovarsSeasonalFuture, 5)
raster25FI = subset(UgandaBiovarsSeasonalFuture, 6)
raster26FI = subset(UgandaBiovarsSeasonalFuture, 7)
raster27FI = subset(UgandaBiovarsSeasonalFuture, 8)
raster28FI = raster("./Uganda08-13/altitude_1KMmedian_MERIT_UgandaClip.tif")
raster29FI = raster("./Uganda08-13/slope_1KMmedian_MERIT_UgandaClip.tif")


#The "I" times 1 puts the raster in active memory
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

raster01F <- raster01FI*1
raster02F <- raster02FI*1
raster03F <- raster03FI*1
raster04F <- raster04FI*1
raster05F <- raster05FI*1
raster06F <- raster06FI*1
raster07F <- raster07FI*1
raster08F <- raster08FI*1
raster09F <- raster09FI*1
raster10F <- raster10FI*1
raster11F <- raster11FI*1
raster12F <- raster12FI*1
raster13F <- raster13FI*1
raster14F <- raster14FI*1
raster15F <- raster15FI*1
raster16F <- raster16FI*1
raster17F <- raster17FI*1
raster18F <- raster18FI*1
raster19F <- raster19FI*1
raster20F <- raster20FI*1
raster21F <- raster21FI*1
raster22F <- raster22FI*1
raster23F <- raster23FI*1
raster24F <- raster24FI*1
raster25F <- raster25FI*1
raster26F <- raster26FI*1
raster27F <- raster27FI*1
raster28F <- raster28FI*1
raster29F <- raster29FI*1

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

proj4string(raster01F) <- crs.geo
proj4string(raster02F) <- crs.geo
proj4string(raster03F) <- crs.geo
proj4string(raster04F) <- crs.geo
proj4string(raster05F) <- crs.geo
proj4string(raster06F) <- crs.geo
proj4string(raster07F) <- crs.geo
proj4string(raster08F) <- crs.geo
proj4string(raster09F) <- crs.geo
proj4string(raster10F) <- crs.geo
proj4string(raster11F) <- crs.geo
proj4string(raster12F) <- crs.geo
proj4string(raster13F) <- crs.geo
proj4string(raster14F) <- crs.geo
proj4string(raster15F) <- crs.geo
proj4string(raster16F) <- crs.geo
proj4string(raster17F) <- crs.geo
proj4string(raster18F) <- crs.geo
proj4string(raster19F) <- crs.geo
proj4string(raster20F) <- crs.geo
proj4string(raster21F) <- crs.geo
proj4string(raster22F) <- crs.geo
proj4string(raster23F) <- crs.geo
proj4string(raster24F) <- crs.geo
proj4string(raster25F) <- crs.geo
proj4string(raster26F) <- crs.geo
proj4string(raster27F) <- crs.geo
proj4string(raster28F) <- crs.geo
proj4string(raster29F) <- crs.geo


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

envF=stack(raster01F,
          raster02F,
          raster03F,
          raster04F,
          raster05F,
          raster06F,
          raster07F,
          raster08F,
          raster09F,
          raster10F,
          raster11F,
          raster12F,
          raster13F,
          raster14F,
          raster15F,
          raster16F,
          raster17F,
          raster18F,
          raster19F,
          raster20F,
          raster21F,
          raster22F,
          raster23F,
          raster24F,
          raster25F,
          raster26F,
          raster27F,
          raster28F,
          raster29F)

names(env) [1]  <- "Bio01_temp_mean_annual"
names(env) [2]  <- "Bio02_temp_mean_diurnal_range"
names(env) [3]  <- "Bio03_isothermality"
names(env) [4]  <- "Bio04_temp_seasonality"
names(env) [5]  <- "Bio05_temp_max_warmestMo"
names(env) [6]  <- "Bio06_temp_min_coldestMo"
names(env) [7]  <- "Bio07_temp_annual_range"
names(env) [8]  <- "Bio08_temp_mean_wettestQ"
names(env) [9]  <- "Bio09_temp_mean_driestQ"  
names(env) [10] <- "Bio10_temp_mean_warmestQ"
names(env) [11] <- "Bio11_temp_mean_coldestQ"
names(env) [12] <- "Bio12_precip_annual"
names(env) [13] <- "Bio13_precip_wettestMo"
names(env) [14] <- "Bio14_precip_driestMo"
names(env) [15] <- "Bio15_precip_seasonality"
names(env) [16] <- "Bio16_precip_wettestQ"
names(env) [17] <- "Bio17_precip_driestQ"
names(env) [18] <- "Bio18_precip_warmestQ"
names(env) [19] <- "Bio19_precip_coldestQ"
names(env) [20] <- "Bio08S_temp_mean_wettestS"
names(env) [21] <- "Bio09S_temp_mean_driestS"  
names(env) [22] <- "Bio10S_temp_mean_warmestS"
names(env) [23] <- "Bio11S_temp_mean_coldestS"
names(env) [24] <- "Bio16S_precip_wettestS"
names(env) [25] <- "Bio17S_precip_driestS"
names(env) [26] <- "Bio18S_precip_warmestS"
names(env) [27] <- "Bio19S_precip_coldestS"
names(env) [28] <- "altitude"
names(env) [29] <- "slope"

#these names have to be the same as the current ones used in building the model.
names(envF) [1]  <- "Bio01_temp_mean_annual"
names(envF) [2]  <- "Bio02_temp_mean_diurnal_range"
names(envF) [3]  <- "Bio03_isothermality"
names(envF) [4]  <- "Bio04_temp_seasonality"
names(envF) [5]  <- "Bio05_temp_max_warmestMo"
names(envF) [6]  <- "Bio06_temp_min_coldestMo"
names(envF) [7]  <- "Bio07_temp_annual_range"
names(envF) [8]  <- "Bio08_temp_mean_wettestQ"
names(envF) [9]  <- "Bio09_temp_mean_driestQ"  
names(envF) [10] <- "Bio10_temp_mean_warmestQ"
names(envF) [11] <- "Bio11_temp_mean_coldestQ"
names(envF) [12] <- "Bio12_precip_annual"
names(envF) [13] <- "Bio13_precip_wettestMo"
names(envF) [14] <- "Bio14_precip_driestMo"
names(envF) [15] <- "Bio15_precip_seasonality"
names(envF) [16] <- "Bio16_precip_wettestQ"
names(envF) [17] <- "Bio17_precip_driestQ"
names(envF) [18] <- "Bio18_precip_warmestQ"
names(envF) [19] <- "Bio19_precip_coldestQ"
names(envF) [20] <- "Bio08S_temp_mean_wettestS"
names(envF) [21] <- "Bio09S_temp_mean_driestS"  
names(envF) [22] <- "Bio10S_temp_mean_warmestS"
names(envF) [23] <- "Bio11S_temp_mean_coldestS"
names(envF) [24] <- "Bio16S_precip_wettestS"
names(envF) [25] <- "Bio17S_precip_driestS"
names(envF) [26] <- "Bio18S_precip_warmestS"
names(envF) [27] <- "Bio19S_precip_coldestS"
names(envF) [28] <- "altitude"
names(envF) [29] <- "slope"

  
print("raster stacks done")

########################################
#Calculate mean of straight lines 
#######################################
StraightMean <- raster::extract(env, spatial.p, fun=mean, na.rm=TRUE)

summary(StraightMean)

StraightMeanDF <- as.data.frame(StraightMean) #DF=data frame

StraightMeanDF$genetdist <- G.table$genetdist #Adds FST into the data frame
StraightMeanDF$geodist <- G.table$geodist #Adds geographic distance into the data frame
head(StraightMeanDF)

#option of trying DPS (another measure of genetic distance)
#StraightMeanDF$DPS <- G.table$DPS

save.image(file = "07a.GeneticRandomForestIterative_prepared.2019-06-18.RData")


########################################
#making initial RF model based on mean of straight lines
#######################################
#with geodist:
Straight_RF1 = randomForest(genetdist ~ geodist +
                            Bio01_temp_mean_annual +
                            Bio02_temp_mean_diurnal_range +
                            Bio03_isothermality +
                            Bio04_temp_seasonality +
                            Bio05_temp_max_warmestMo +
                            Bio06_temp_min_coldestMo +
                            Bio07_temp_annual_range +
                            Bio08_temp_mean_wettestQ +
                            Bio09_temp_mean_driestQ +  
                            Bio10_temp_mean_warmestQ +
                            Bio11_temp_mean_coldestQ +
                            Bio12_precip_annual +
                            Bio13_precip_wettestMo +
                            Bio14_precip_driestMo +
                            Bio15_precip_seasonality +
                            Bio16_precip_wettestQ +
                            Bio17_precip_driestQ +
                            Bio18_precip_warmestQ +
                            Bio19_precip_coldestQ +
                            Bio08S_temp_mean_wettestS +
                            Bio09S_temp_mean_driestS +  
                            Bio10S_temp_mean_warmestS +
                            Bio11S_temp_mean_coldestS +
                            Bio16S_precip_wettestS +
                            Bio17S_precip_driestS +
                            Bio18S_precip_warmestS +
                            Bio19S_precip_coldestS +
                            altitude +
                            slope, importance=TRUE, na.action=na.omit, data=StraightMeanDF)

Straight_RF1 #gives you the % variance explained
varImpPlot(Straight_RF1) #gives a plot of which variables are most important

StraightPred1 <- predict(env, Straight_RF1) #doesn't work because the model has an extra variable, geodist that the raster, env, does not
#plot(StraightPred1)
#pred.cond1 <- 1/StraightPred1 #build conductance surface
#plot(pred.cond1)

#without geodist:
Straight_RF = randomForest(genetdist ~
                             Bio01_temp_mean_annual +
                             Bio02_temp_mean_diurnal_range +
                             Bio03_isothermality +
                             Bio04_temp_seasonality +
                             Bio05_temp_max_warmestMo +
                             Bio06_temp_min_coldestMo +
                             Bio07_temp_annual_range +
                             Bio08_temp_mean_wettestQ +
                             Bio09_temp_mean_driestQ +  
                             Bio10_temp_mean_warmestQ +
                             Bio11_temp_mean_coldestQ +
                             Bio12_precip_annual +
                             Bio13_precip_wettestMo +
                             Bio14_precip_driestMo +
                             Bio15_precip_seasonality +
                             Bio16_precip_wettestQ +
                             Bio17_precip_driestQ +
                             Bio18_precip_warmestQ +
                             Bio19_precip_coldestQ +
                             Bio08S_temp_mean_wettestS +
                             Bio09S_temp_mean_driestS +  
                             Bio10S_temp_mean_warmestS +
                             Bio11S_temp_mean_coldestS +
                             Bio16S_precip_wettestS +
                             Bio17S_precip_driestS +
                             Bio18S_precip_warmestS +
                             Bio19S_precip_coldestS +
                             altitude +
                             slope, importance=TRUE, na.action=na.omit, data=StraightMeanDF)

Straight_RF #gives you the % variance explained
varImpPlot(Straight_RF) #gives a plot of which variables are most important


print("first prediction resistance surface done")

#Save both results of variable imnportance with and without geodist:
dev.off()
pdf(file="VariableImportance_08-13.pdf",height=6,width=10)
par(mfrow=c(1,2))
varImpPlot(Straight_RF, main="Random Forest variance importance (without geo. dist.)")
varImpPlot(Straight_RF1, main="Random Forest variance importance (with geo. dist.)")
dev.off()


###########################
#Make current prediction:
###########################
StraightPred <- predict(env, Straight_RF) 
plot(StraightPred,frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)

pred.cond <- 1/StraightPred #build conductance surface
plot(pred.cond,frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)

###########################
#Make future prediction:
###########################
StraightPredF <- predict(envF, Straight_RF) 
plot(StraightPredF,frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)

pred.condF <- 1/StraightPredF #build conductance surface
plot(pred.condF,frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)

###########################
# to save a raster file, for example:
###########################
writeRaster(StraightPred, "RF_ResistanceSurface_2008-13.asc", format = "ascii", overwrite=T)
writeRaster(StraightPredF, "RF_ResistanceSurface_RCP45.asc", format = "ascii", overwrite=T)
writeRaster(pred.cond, "RF_ConnectivitySurface_2008-13.asc", format = "ascii", overwrite=T)
writeRaster(pred.condF, "RF_ConnectivitySurface_RCP45.asc", format = "ascii", overwrite=T)


###############################################
#MAP some of the rasters to visualize
###############################################
#load("~/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-RF-Genetic-Model/07a.GeneticRandomForestIterative.2019-06-18.RData")

#Habitat suitability:
FAO <- raster("~/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/fuscipesf/fuscipesf/dblbnd.adf") #FAO map
proj4string(FAO) <- crs.geo 
plot(FAO)

#water shape of Uganda:
waterFull <- shapefile("~/Dropbox/Caccone_Aksoy/QGIS/from\ www.diva-gis.org/waterbodies_africa/waterbodies_africa.dbf")
spWater <- spTransform(waterFull,crs.geo)
waterFull <- spWater
plot(waterFull, col="black", legend=FALSE)
waterUg <- crop(waterFull,extent(28.6, 35.4, -1.5, 4.73))
plot(waterUg, col="black", legend=FALSE)

#countries of Africa:
boundriesFull <- shapefile("~/Dropbox/Caccone_Aksoy/QGIS/from\ www.diva-gis.org/AfricanCountries/AfricanCountires.shp")
spBoundaries <- spTransform(boundriesFull,crs.geo)
boundriesFull<- spBoundaries
plot(boundriesFull, legend=FALSE)
boundariesUg <- crop(boundriesFull,extent(28.6, 35.4, -1.5, 4.73))
plot(boundariesUg, legend=FALSE)

colfunc1a <- colorRampPalette(c("white","light yellow","orange","green", "blue","dark blue","black"))
colfunc1b <- colorRampPalette(c("white","light yellow","orange","green"))
colfunc1c <- colorRampPalette(c("light yellow","orange","green","blue"))
colfunc2 <- colorRampPalette(c("dark red","gold","light yellow","white"))
colfunc3 <- colorRampPalette(c("white","dark green"))
colfunc4 <- colorRampPalette(c("light green","white","light yellow", "gold","dark red","red"))
colfunc5 <- colorRampPalette(c("light green","white","gold","dark red","red"))
colfunc6 <- colorRampPalette(c("dark red","gold","white","light green","blue"))
colfunc7 <- colorRampPalette(c("blue","light green","white","gold","dark red"))



###########################
#plot current and future predictions side by side:
###########################
par(mfrow=c(2,1))
plot(StraightPred,frame.plot=F,axes=F,box=F,add=F,legend=T,col=colfunc1a(20))
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)
plot(StraightPredF,frame.plot=F,axes=F,box=F,add=F,legend=T,col=colfunc1b(20))
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)

###########################
#plot current and future connectivity (inverse of predicted Fst) side by side:
###########################
par(mfrow=c(2,1))
plot(pred.cond,frame.plot=F,axes=F,box=F,add=F,legend=T,col=colfunc1a(20),main="Predicted Genetic Connectivity 2008-2013")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)
plot(pred.condF,frame.plot=F,axes=F,box=F,add=F,legend=T,col=colfunc1c(20),main="Predicted Genetic Connectivity 2041-2060 (RCP 4.5)")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)


###########################
#Plot the difference between the current and future:
###########################
rf.c<-pred.cond
rf.f<-pred.condF
proj4string(rf.c)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(rf.f)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#Use overlay to determine diff between future and current
rf.dif<-overlay(rf.f,rf.c,fun=function(r1, r2){return(r1-r2)})

#Plot with simple command
plot(rf.dif,frame.plot=F,axes=F,box=F,add=F,legend=T,col=colfunc6(16),main="Predicted Change in Genetic Connectivity")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)

###########################
#save project:
###########################
save.image(file = "07a.GeneticRandomForest.2019-06-22.RData")


###########################
#plot current and future top ranking variable:
###########################

### 1st ranking variable:
#pdf(file="Bio03_isothermality_rank1.pdf")
par(mfrow=c(2,1))
plot(env$Bio03_isothermality,main="Rank 1: Bio03 isothermality 2008-2013",col=colfunc1(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
map(interior=T,xlim=c(28.6, 35.4),ylim=c(-1.5, 4.73),add=T)
plot(waterUg,add=T,col="black",legend=F)
plot(envF$Bio03_isothermality,main="Rank 1: Bio03 isothermality RP 4.5 forecast",col=colfunc1(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
map(interior=T,add=T)
plot(waterUg,add=T,col="black",legend=F)
#dev.off()
#Calculate and map the difference
is.c<-env$Bio03_isothermality
is.f<-envF$Bio03_isothermality
proj4string(is.c)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(is.f)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#Use overlay to determine diff between future and current
is.dif<-overlay(is.f,is.c,fun=function(r1, r2){return(r1-r2)})
#Plot with simple command
plot(is.dif,col=colfunc4(10),frame.plot=F,axes=F,box=F,add=F,legend=T,main="Change in Isothermality")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)

### 2nd ranking variable:
pdf(file="Bio15_precip_seasonality_rank2.pdf")
par(mfrow=c(2,1))
plot(env$Bio15_precip_seasonality,main="Rank 2: Bio15 precip seasonality 2008-2013",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
map(interior=T,region="Uganda",add=T)
plot(waterUg,add=T,col="black",legend=F)
plot(envF$Bio15_precip_seasonality,main="Rank 2: Bio15 precip seasonality RP 4.5 forecast",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
map(interior=T,add=T)
plot(waterUg,add=T,col="black",legend=F)
dev.off()
#Calculate and map the difference
prp.c<-env$Bio15_precip_seasonality
prp.f<-envF$Bio15_precip_seasonality
proj4string(prp.c)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(prp.f)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#Use overlay to determine diff between future and current
prp.dif<-overlay(prp.f,prp.c,fun=function(r1, r2){return(r1-r2)})
#Plot with simple command
plot(prp.dif,col=colfunc2(16),frame.plot=F,axes=F,box=F,add=F,legend=T,main="Change in Precipitation Seasonality")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)

### 3rd ranking variable:
#pdf(file="Bio08S_temp_mean_wettestS_rank3.pdf")
par(mfrow=c(2,1))
plot(env$Bio08S_temp_mean_wettestS,main="Rank 3: Bio08S temp mean wettest S 2008-2013",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(boundariesUg,add=T)
plot(waterUg,add=T,col="black",legend=F)
plot(envF$Bio08S_temp_mean_wettestS,main="Rank 3: Bio08S temp mean wettest S RP 4.5 forecast",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(boundariesUg,add=T)
plot(waterUg,add=T,col="black",legend=F)
#dev.off()
#Calculate and map the difference
stemp.c<-env$Bio08S_temp_mean_wettestS
stemp.f<-envF$Bio08S_temp_mean_wettestS
proj4string(stemp.c)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(stemp.f)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#Use overlay to determine diff between future and current
stemp.dif<-overlay(stemp.f,stemp.c,fun=function(r1, r2){return(r1-r2)})
#Plot with simple command
plot(stemp.dif,col=colfunc4(10),frame.plot=F,axes=F,box=F,add=F,legend=T,main="Change in Mean Temp. in Wettest Season")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)

### 4th ranking variable:
#pdf(file="Bio14_precip_driestMo_rank4.pdf")
par(mfrow=c(2,1))
plot(env$Bio14_precip_driestMo,main="Rank 4: Bio14 precip driest Mo 2008-2013",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(boundariesUg,add=T)
plot(waterUg,add=T,col="black",legend=F)
plot(envF$Bio14_precip_driestMo,main="Rank 4: Bio14 precip driest Mo RP 4.5 forecast",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
plot(boundariesUg,add=T)
plot(waterUg,add=T,col="black",legend=F)
#dev.off()
#Calculate and map the difference
stemp.c<-env$Bio14_precip_driestMo
stemp.f<-envF$Bio14_precip_driestMo
proj4string(stemp.c)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
proj4string(stemp.f)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
#Use overlay to determine diff between future and current
prcp.dif<-overlay(stemp.f,stemp.c,fun=function(r1, r2){return(r1-r2)})
#Plot with simple command
plot(prcp.dif,col=colfunc5(10),frame.plot=F,axes=F,box=F,add=F,legend=T,main="Change in Precipitation in Driest Month")
plot(waterUg,add=T,col="black",legend=F)
plot(boundariesUg,add=T)


#add lines (data units) with FST and mean under straight lines used in the model
#pdf(file="StraightLines.pdf")
plot(env$altitude,main="pairwise distance lines",col=colfunc2(25),frame.plot=F,axes=F,box=F,add=F,legend=F)
map(interior=T,add=T)
plot(waterUg,add=T,col="black",legend=F)
plot(p,add=T,col="blue")
#dev.off()



###########################
#save project:
###########################
save.image(file = "07a.GeneticRandomForest.2019-06-22.RData")


###########################
# Validate with a subsampling approach
###########################

NumPairs = nrow(StraightMeanDF)
Training = NumPairs * 0.7
TrainingInt = round(Training)
TrainingPairs = sample(1:NumPairs, TrainingInt, replace = FALSE)
StraightMeanDF.train = StraightMeanDF[TrainingPairs,]
StraightMeanDF.valid = StraightMeanDF[-TrainingPairs,]

Straight_RF2 = randomForest(genetdist ~
                              Bio01_temp_mean_annual +
                              Bio02_temp_mean_diurnal_range +
                              Bio03_isothermality +
                              Bio04_temp_seasonality +
                              Bio05_temp_max_warmestMo +
                              Bio06_temp_min_coldestMo +
                              Bio07_temp_annual_range +
                              Bio08_temp_mean_wettestQ +
                              Bio09_temp_mean_driestQ +  
                              Bio10_temp_mean_warmestQ +
                              Bio11_temp_mean_coldestQ +
                              Bio12_precip_annual +
                              Bio13_precip_wettestMo +
                              Bio14_precip_driestMo +
                              Bio15_precip_seasonality +
                              Bio16_precip_wettestQ +
                              Bio17_precip_driestQ +
                              Bio18_precip_warmestQ +
                              Bio19_precip_coldestQ +
                              Bio08S_temp_mean_wettestS +
                              Bio09S_temp_mean_driestS +  
                              Bio10S_temp_mean_warmestS +
                              Bio11S_temp_mean_coldestS +
                              Bio16S_precip_wettestS +
                              Bio17S_precip_driestS +
                              Bio18S_precip_warmestS +
                              Bio19S_precip_coldestS +
                              altitude +
                              slope, importance=TRUE, na.action=na.omit, data=StraightMeanDF.train)


#why aren't these the same value?
Straight_rsq = tail(Straight_RF2$rsq ,1 ) 
Straight_rsq_fromEquation = 1 - (tail(Straight_RF2$mse , 1) / var(StraightMeanDF.train$genetdist))     
write.table(Straight_rsq, "RSQ_Table.txt")
write.table(Straight_rsq_fromEquation, "RSQ_fromEquation_Table.txt")

Straight_mse = tail(Straight_RF2$mse ,1 )  
Straight_mse_fromEquation = mean ((predict(Straight_RF2, StraightMeanDF.train) - StraightMeanDF.train$genetdist)^2) 
Straight_mse2 = mean ((predict(Straight_RF2, StraightMeanDF.valid) - StraightMeanDF.valid$genetdist)^2)
write.table(Straight_mse, "MSE_Table.txt")
write.table(Straight_mse_fromEquation, "MSE_fromEquation_Table.txt")
write.table(Straight_mse2, "MSE2_Table.txt")

#save these as variables and save to a table
cor1 = cor(Straight_RF2$predict, StraightMeanDF.train$genetdist)
cor2 = cor(predict(Straight_RF2, StraightMeanDF.valid), StraightMeanDF.valid$genetdist)
write.table(cor1, "InternalValidation.txt")
write.table(cor2, "ExternalValidation.txt")



#Plot observed vs predicted
Tdata <- StraightMeanDF.train$genetdist
Vdata <- StraightMeanDF.valid$genetdist
Idata <- StraightMeanDF$genetdist
PredictedIdata <- Straight_RF$predicted
PredictedTdata <-raster::predict(Straight_RF2, StraightMeanDF.train)
PredictedVdata <- raster::predict(Straight_RF2, StraightMeanDF.valid)

#set up PDF file:
pdf("Observed_vs_Predicted_withoutGeoDist.pdf", 4, 10,paper="letter")
par(mfrow=c(3,1))

#Internal validation:
linearModI <- lm(Idata ~ PredictedIdata)
cor(Idata,PredictedIdata) #0.8567832
summary(linearModI) #Residual standard error: 0.04808 on 1643 degrees of freedom, Multiple R-squared:  0.7341,Adjusted R-squared:  0.7339 , F-statistic:  4535 on 1 and 1643 DF,  p-value: < 2.2e-16
plot(Idata, PredictedIdata, xlab="Training Data",ylab="Predicted from Out-of-Bag Samples",main="Internal Validation")
abline(lm(Idata ~ PredictedIdata), col="red")
text(0.52, 0.22, "R-squared:  0.7341 \n p-value: < 2.2e-16", col = "red")

#linear model of Training data:
linearModT <- lm(Tdata ~ PredictedTdata)
cor(Tdata,PredictedTdata) #0.9766753
summary(linearModT) #Residual standard error: 0.01982 on 1150 degrees of freedom, Multiple R-squared:  0.9539,Adjusted R-squared:  0.9539, F-statistic: 2.379e+04 on 1 and 1150 DF,  p-value: < 2.2e-16
plot(Tdata, PredictedTdata, xlab="Training Data",ylab="Predicted Training Data",main="Validation based on 70% Training Data")
abline(lm(Tdata ~ PredictedTdata), col="red")
text(0.52, 0.22, "R-squared:  0.9539 \n p-value: < 2.2e-16", col = "red")

#linear model of Validation data:
linearModV <- lm(Vdata ~ PredictedVdata)
cor(Vdata,PredictedVdata) #0.8501442
summary(linearModV) #Residual standard error: 0.05031 on 491 degrees of freedom, Multiple R-squared:  0.7227,Adjusted R-squared:  0.7222 , F-statistic:  1280 on 1 and 491 DF,  p-value: < 2.2e-16
plot(Vdata, PredictedVdata, xlab="Validation Data",ylab="Predicted Validation Data",main="Validation based on 30% Validation Data")
abline(lm(Vdata ~ PredictedVdata), col="red")
text(0.52, 0.22, "R-squared:  0.7227 \n p-value: < 2.2e-16", col = "red")


dev.off()

StraightPred2 <- predict(env, Straight_RF2) #if i want to see the  model based on 70% of the data without geo. dist.


###########################
# Validate with a subsampling approach including geographic distance:
###########################

NumPairs = nrow(StraightMeanDF)
Training = NumPairs * 0.7
TrainingInt = round(Training)
TrainingPairs = sample(1:NumPairs, TrainingInt, replace = FALSE)
StraightMeanDF.train = StraightMeanDF[TrainingPairs,]
StraightMeanDF.valid = StraightMeanDF[-TrainingPairs,]

Straight_RF3 = randomForest(genetdist ~ geodist +
                              Bio01_temp_mean_annual +
                              Bio02_temp_mean_diurnal_range +
                              Bio03_isothermality +
                              Bio04_temp_seasonality +
                              Bio05_temp_max_warmestMo +
                              Bio06_temp_min_coldestMo +
                              Bio07_temp_annual_range +
                              Bio08_temp_mean_wettestQ +
                              Bio09_temp_mean_driestQ +  
                              Bio10_temp_mean_warmestQ +
                              Bio11_temp_mean_coldestQ +
                              Bio12_precip_annual +
                              Bio13_precip_wettestMo +
                              Bio14_precip_driestMo +
                              Bio15_precip_seasonality +
                              Bio16_precip_wettestQ +
                              Bio17_precip_driestQ +
                              Bio18_precip_warmestQ +
                              Bio19_precip_coldestQ +
                              Bio08S_temp_mean_wettestS +
                              Bio09S_temp_mean_driestS +  
                              Bio10S_temp_mean_warmestS +
                              Bio11S_temp_mean_coldestS +
                              Bio16S_precip_wettestS +
                              Bio17S_precip_driestS +
                              Bio18S_precip_warmestS +
                              Bio19S_precip_coldestS +
                              altitude +
                              slope, importance=TRUE, na.action=na.omit, data=StraightMeanDF.train)


#why aren't these the same value?
Straight_rsq = tail(Straight_RF3$rsq ,1 ) 
Straight_rsq_fromEquation = 1 - (tail(Straight_RF3$mse , 1) / var(StraightMeanDF.train$genetdist))     
write.table(Straight_rsq, "RSQ_Table.txt")
write.table(Straight_rsq_fromEquation, "RSQ_fromEquation_Table.txt")

Straight_mse = tail(Straight_RF3$mse ,1 )  
Straight_mse_fromEquation = mean ((predict(Straight_RF3, StraightMeanDF.train) - StraightMeanDF.train$genetdist)^2) 
Straight_mse2 = mean ((predict(Straight_RF3, StraightMeanDF.valid) - StraightMeanDF.valid$genetdist)^2)
write.table(Straight_mse, "MSE_Table.txt")
write.table(Straight_mse_fromEquation, "MSE_fromEquation_Table.txt")
write.table(Straight_mse2, "MSE2_Table.txt")

#save these as variables and save to a table
cor1 = cor(Straight_RF3$predict, StraightMeanDF.train$genetdist)
cor2 = cor(predict(Straight_RF3, StraightMeanDF.valid), StraightMeanDF.valid$genetdist)
write.table(cor1, "InternalValidation.txt")
write.table(cor2, "ExternalValidation.txt")



#Plot observed vs predicted
Tdata <- StraightMeanDF.train$genetdist
Vdata <- StraightMeanDF.valid$genetdist
Idata <- StraightMeanDF$genetdist
PredictedIdata <- Straight_RF1$predicted
PredictedTdata <-raster::predict(Straight_RF3, StraightMeanDF.train)
PredictedVdata <- raster::predict(Straight_RF3, StraightMeanDF.valid)

#Start PDF file:
pdf("Observed_vs_Predicted_withGeoDist.pdf", 4, 10,paper="letter")
par(mfrow=c(3,1))

#Internal validation:
linearModI <- lm(Idata ~ PredictedIdata)
cor(Idata,PredictedIdata) #0.8964952
summary(linearModI) #Residual standard error: 0.04131 on 1643 degrees of freedom, Multiple R-squared:  0.8037,Adjusted R-squared:  0.8036, F-statistic:  6727 on 1 and 1643 DF,  p-value: < 2.2e-16
plot(Idata, PredictedIdata, xlab="Training Data",ylab="Predicted from Out-of-Bag Samples",main="Internal Validation with Geo. Dist.")
abline(lm(Idata ~ PredictedIdata), col="red")
text(0.52, 0.22, "R-squared:  0.8037 \np-value: < 2.2e-16", col = "red")

#linear model of Training data:
linearModT <- lm(Tdata ~ PredictedTdata)
cor(Tdata,PredictedTdata) #0.9840072
summary(linearModT) #Residual standard error: 0.01706 on 1150 degrees of freedom, Multiple R-squared:  0.9683,Adjusted R-squared:  0.9682 , F-statistic: 3.509e+04 on 1 and 1150 DF,  p-value: < 2.2e-16
plot(Tdata, PredictedTdata, xlab="Training Data",ylab="Predicted Training Data",main="Validation based on 70% Training Data with Geo. Dist.")
abline(lm(Tdata ~ PredictedTdata), col="red")
text(0.52, 0.22, "R-squared:  0.9683 \np-value: < 2.2e-16", col = "red")

#linear model of Validation data:
linearModV <- lm(Vdata ~ PredictedVdata)
cor(Vdata,PredictedVdata) #0.8629154
summary(linearModV) #Residual standard error: 0.04402 on 491 degrees of freedom, Multiple R-squared:  0.7446,Adjusted R-squared:  0.7441 , F-statistic:  1432 on 1 and 491 DF,  p-value: < 2.2e-16
plot(Vdata, PredictedVdata, xlab="Validation Data",ylab="Predicted Validation Data",main="Validation based on 30% Validation Data with Geo. Dist.")
abline(lm(Vdata ~ PredictedVdata), col="red")
text(0.52, 0.22, "R-squared:  0.7446 \np-value: < 2.2e-16", col = "red")

dev.off()



###########################
#save project:
###########################
save.image(file = "07a.GeneticRandomForest.2019-06-22.RData")


#write this as tif and remove it

print("first prediction resistance surface done")



