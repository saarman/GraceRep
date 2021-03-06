
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

##################################################################################
#Go back to look at Random Forest results:
##################################################################################
setwd("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-RF-Genetic-Model")
load("~/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-RF-Genetic-Model/03a.GeneticRandomForest.2019-04-06.RData")
library("randomForest")

#The prediction, of resistance, and of connectivity:
plot(StraightPred)
StraightPredConnect <- 1/StraightPred #build conductance surface
plot(StraightPredConnect)

#look at both results with and without geodist:
dev.off()
pdf(file="PCA-ddRAD-v15-April3.pdf",height=6,width=10)
par(mfrow=c(1,2))
varImpPlot(Straight_RF, main="Random Forest variance importance (without geo. dist.)")
varImpPlot(Straight_RF1, main="Random Forest variance importance (with geo. dist.)")
dev.off()

##################################################################################
#Make bivariate.map
##################################################################################
scaledColor <- c("#AFCD93", "#FAF080", "#FACE7F", "#EEA32D", "#DF512A", "#D1232A", "#822024")
waterColor <- "#9AD5E6"

library(classInt)
library(raster)
library(rgdal)
library(dismo)
library(XML)
library(maps)
library(sp)
library(classInt)

####################################################
#import files, include FAO and maxHS raster file
####################################################

PRED <- raster("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-spatial/Gff-Uganda-RF-Genetic-Model/RF_ConnectivitySurfaceScaled.asc")
proj4string(PRED) <- crs.geo 
plot(PRED)

FAO <- raster("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/fuscipesf/fuscipesf/dblbnd.adf") #FAO map
proj4string(FAO) <- crs.geo 
plot(FAO)

maryHS <- raster("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/HS-maps/PSN_raw/Glossina_fuscipes_fuscipes.asc") 
proj4string(maryHS) <- crs.geo 
plot(maryHS)

meanHS <- raster("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/HS-maps/PSN_FAO_average/psnfao_mean/dblbnd.adf")
proj4string(meanHS) <- crs.geo 
plot(meanHS)

maxHS <- raster("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/HS-maps/PSN_FAO_max/psnfao_max/dblbnd.adf")  #max of FAO and updated habitat suitability model
proj4string(maxHS) <- crs.geo 
plot(maxHS)

blendHS <- raster("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/HS-maps/PSN_FAO_blend/psnfao_blend/dblbnd.adf") #max of FAO and updated habitat suitability
proj4string(blendHS) <- crs.geo 
plot(blendHS)

waterFull <- shapefile("/Users/nsaarman1/Dropbox/Caccone_Aksoy/QGIS/from\ www.diva-gis.org/waterbodies_africa/waterbodies_africa.dbf") 
spWater <- spTransform(waterFull,crs.geo)
waterFull <- spWater
plot(waterFull, col="black", legend=FALSE)

Uganda <- shapefile("/Users/nsaarman1/Dropbox/Caccone_Aksoy/Glossina-Mary/R-tools-patches/UGA_adm/UGA_adm0.dbf")
spUganda <- spTransform(Uganda,crs.geo)
Uganda <- spUganda
plot(Uganda)

####################################################
# Compare different maps:
####################################################
dev.off()
plot(FAO,main="FAO habitat suitability map from 2001",frame.plot=F,axes=F,box=F,add=F,legend=F)
map(interior=T,add=T)
plot(spWater,add=T,col="black",legend=F)

plot(blendHS,main="Updated FAO habitat suitability map from 2018",frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)
plot(spWater,add=T,col="black",legend=F)

plot(PRED,main="Random Forest Predicted Connectivity",frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)
plot(spWater,add=T,col="black",legend=F)

####################################
#crop each HS to Uganda
####################################
FAOcrop <- crop(FAO,Uganda)
#rescale to 0 to 1
r <- FAOcrop
r.min = cellStats(r, "min")
r.max = cellStats(r, "max")
r.scale <- ((r - r.min) / (r.max - r.min))
FAOcrop <- r.scale

######################
blendHScrop <- crop(blendHS,Uganda)
#rescale to 0 to 1
r <- blendHScrop
r.min = cellStats(r, "min")
r.max = cellStats(r, "max")
r.scale <- ((r - r.min) / (r.max - r.min))
blendHScrop <- r.scale
######################

PREDcrop <- crop(PRED,Uganda)
#rescale to 0 to 1
r <- PREDcrop
r.min = cellStats(r, "min")
r.max = cellStats(r, "max")
r.scale <- ((r - r.min) / (r.max - r.min))
PREDcrop <- r.scale
##########


#############################################
# Compare different maps after cropping
dev.off()
plot(FAOcrop,main="FAO habitat suitability map from 2001",frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)
plot(spWater,add=T,col="black",legend=F)

plot(blendHScrop,main="FAO habitat suitability map updated 2018",frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)
plot(spWater,add=T,col="black",legend=F)

plot(PREDcrop,main="Random Forest Predicted Connectivity",frame.plot=F,axes=F,box=F,add=F,legend=T)
map(interior=T,add=T)
plot(spWater,add=T,col="black",legend=F)

#############################################
#create a color matrix for Bivariate map:
#############################################
colmat<-function(nquantiles=10, upperleft=rgb(0,150,235, maxColorValue=255), upperright=rgb(130,0,80, maxColorValue=255), bottomleft="grey", bottomright=rgb(255,230,15, maxColorValue=255), xlab="x label", ylab="y label"){
  my.data<-seq(0,1,.01)
  my.class<-classIntervals(my.data,n=nquantiles,style="quantile")
  my.pal.1<-findColours(my.class,c(upperleft,bottomleft))
  my.pal.2<-findColours(my.class,c(upperright, bottomright))
  col.matrix<-matrix(nrow = 101, ncol = 101, NA)
  for(i in 1:101){
    my.col<-c(paste(my.pal.1[i]),paste(my.pal.2[i]))
    col.matrix[102-i,]<-findColours(my.class,my.col)}
  plot(c(1,1),pch=19,col=my.pal.1, cex=0.5,xlim=c(0,1),ylim=c(0,1),frame.plot=F, xlab=xlab, ylab=ylab,cex.lab=1.3)
  for(i in 1:101){
    col.temp<-col.matrix[i-1,]
    points(my.data,rep((i-1)/100,101),pch=15,col=col.temp, cex=1)}
  seqs<-seq(0,100,(100/nquantiles))
  seqs[1]<-1
  col.matrix<-col.matrix[c(seqs), c(seqs)]}

col.matrix<-colmat(nquantiles=10, xlab="genetic connectivity", ylab="habitat suitability")

#############################################
#Build bivariate map function:
#############################################
bivariate.map<-function(rasterx, rastery, colormatrix=col.matrix, nquantiles=10){
  quanmean<-getValues(rasterx)
  temp<-data.frame(quanmean, quantile=rep(NA, length(quanmean)))
  brks<-with(temp, quantile(temp,na.rm=T, probs = c(seq(0,1,1/nquantiles))))
  r1<-within(temp, quantile <- cut(quanmean, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr<-data.frame(r1[,2]) 
  quanvar<-getValues(rastery)
  temp<-data.frame(quanvar, quantile=rep(NA, length(quanvar)))
  brks<-with(temp, quantile(temp,na.rm=TRUE, probs = c(seq(0,1,1/nquantiles))))
  r2<-within(temp, quantile <- cut(quanvar, breaks = brks, labels = 2:length(brks),include.lowest = TRUE))
  quantr2<-data.frame(r2[,2])
  as.numeric.factor<-function(x) {as.numeric(levels(x))[x]}
  col.matrix2<-colormatrix
  cn<-unique(colormatrix)
  for(i in 1:length(col.matrix2)){
    ifelse(is.na(col.matrix2[i]),col.matrix2[i]<-1,col.matrix2[i]<-which(col.matrix2[i]==cn)[1])}
  cols<-numeric(length(quantr[,1]))
  for(i in 1:length(quantr[,1])){
    a<-as.numeric.factor(quantr[i,1])
    b<-as.numeric.factor(quantr2[i,1])
    cols[i]<-as.numeric(col.matrix2[b,a])}
  r<-rasterx
  r[1:length(r)]<-cols
  return(r)}



#############################################
# Use the bivariate.map function:
#############################################

bivmap<-bivariate.map(PREDcrop,FAOcrop,colormatrix=col.matrix, nquantiles=10)
# Note: go back to the function that generates color matrices and try out new color schemes.

# Plot the bivariate map:
plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
map(interior=T,add=T)

#plot all three:
dev.off()
par(mfrow=c(3,1))
plot(PREDcrop,main="Random Forest Predicted Connectivity") #, col=scaledColor, breaks=c(0,.1,.25,.5,.75,.9,.95,1))
map(interior=T,add=T)
#plot(Uganda, add=T)
plot(waterFull,add=T,col="black",legend=F)

plot(maxHScrop,main="Improved habitat suitability map") #, col=scaledColor, breaks=c(0,.1,.25,.5,.75,.9,.95,1))
map(interior=T,add=T)
#plot(Uganda, add=T)
plot(waterFull,add=T,col="black",legend=F)

plot(bivmap,main="Genetic Connectivity vs Habitat Suitability",legend=F, col=as.vector(col.matrix))
map(interior=T,add=T)
plot(waterFull,add=T,col="black",legend=F)


########################################################################################################
#Now try the bivariate map with maxHS layer:
# Use the bivariate.map function:
bivmap<-bivariate.map(PREDcrop,maxHScrop,colormatrix=col.matrix, nquantiles=10)
# Note: go back to the function that generates color matrices and try out new color schemes.

# Plot the bivariate map:
plot(bivmap,frame.plot=F,axes=F,box=F,add=F,legend=F,col=as.vector(col.matrix))
map(interior=T,add=T)

#plot all three:
dev.off()
par(mfrow=c(3,1))
plot(PREDcrop,main="Random Forest Predicted Connectivity") #, col=scaledColor, breaks=c(0,.1,.25,.5,.75,.9,.95,1))
map(interior=T,add=T)
#plot(Uganda, add=T)
plot(waterFull,add=T,col="black",legend=F)

plot(maxHScrop,main="Improved habitat suitability map") #, col=scaledColor, breaks=c(0,.1,.25,.5,.75,.9,.95,1))
map(interior=T,add=T)
#plot(Uganda, add=T)
plot(waterFull,add=T,col="black",legend=F)

plot(bivmap,main="Genetic Connectivity vs Habitat Suitability",legend=F, col=as.vector(col.matrix))
map(interior=T,add=T)
plot(waterFull,add=T,col="black",legend=F)
