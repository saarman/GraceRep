
#first thing you do is load the library that is for this analysis (as in adegenet package)
library("adegenet")
library("hierfstat")
library("reshape2")
library("plotly")
library("dplyr")
library("classInt")
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
library("SDraw")

#set working directory
setwd("~/Dropbox/Caccone_Aksoy/Glossina-spatial/Gpd-RF-Genetic-Model")

#set coordinate system for mapping
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")



#load and edit genind object
##############################################
#load genepop file as a genind object and population info file
################################
myData <- read.genepop("input_Gpd_KenTza_11loci_700indv_NPS.gen", ncode = 3)

pop_info <- read.table("input_Gpd_KenTza_11loci_700indv_NPS_info_corrected.csv",header=T,sep=",")

#add the name of the populations to the geneind object
myData@pop <- pop_info$siteID
myData@other$a_xy <- pop_info[,c("approxLong","approxLat")]
myData@other$e_xy <- pop_info[,c("corLong","corLat")]

##################################
#define regions and clusters
##################################
pop_info$region<- as.character(pop_info$siteID)
pop_info$region[pop_info$siteID=="01_KAP"] <- "northwest"
pop_info$region[pop_info$siteID=="02_OTT"] <- "northwest"
pop_info$region[pop_info$siteID=="03_RUM"] <- "northwest"
pop_info$region[pop_info$siteID=="04_GVR"] <- "southwest"
pop_info$region[pop_info$siteID=="05_MRB"] <- "southwest"
pop_info$region[pop_info$siteID=="06_MRT"] <- "southwest"
pop_info$region[pop_info$siteID=="07_NBS"] <- "southwest"
pop_info$region[pop_info$siteID=="08_FGT"] <- "southwest"
pop_info$region[pop_info$siteID=="09_IKR"] <- "southwest"
pop_info$region[pop_info$siteID=="10_GTR"] <- "southwest"
pop_info$region[pop_info$siteID=="11_KLM"] <- "southwest"
pop_info$region[pop_info$siteID=="12_MSN"] <- "southwest"
pop_info$region[pop_info$siteID=="13_MSS"] <- "southwest"
pop_info$region[pop_info$siteID=="14_NGK"] <- "southwest"
pop_info$region[pop_info$siteID=="15_NGU"] <- "central"
pop_info$region[pop_info$siteID=="16_MNP"] <- "east"
pop_info$region[pop_info$siteID=="17_HND"] <- "east"
pop_info$region[pop_info$siteID=="18_AMR"] <- "east"
pop_info$region[pop_info$siteID=="19_CNP"] <- "east"
pop_info$region[pop_info$siteID=="20_KIB"] <- "east"
pop_info$region[pop_info$siteID=="21_TSW"] <- "east"
pop_info$region[pop_info$siteID=="22_KIN"] <- "east"
pop_info$region[pop_info$siteID=="23_SHI"] <- "east"
pop_info$region[pop_info$siteID=="24_SHT"] <- "east"
myData@other$region <- pop_info$region

pop_info$cluster <- as.character(pop_info$siteID)
pop_info$cluster[pop_info$siteID=="01_KAP"] <- "west"
pop_info$cluster[pop_info$siteID=="02_OTT"] <- "west"
pop_info$cluster[pop_info$siteID=="03_RUM"] <- "west"
pop_info$cluster[pop_info$siteID=="04_GVR"] <- "west"
pop_info$cluster[pop_info$siteID=="05_MRB"] <- "west"
pop_info$cluster[pop_info$siteID=="06_MRT"] <- "west"
pop_info$cluster[pop_info$siteID=="07_NBS"] <- "west"
pop_info$cluster[pop_info$siteID=="08_FGT"] <- "west"
pop_info$cluster[pop_info$siteID=="09_IKR"] <- "west"
pop_info$cluster[pop_info$siteID=="10_GTR"] <- "west"
pop_info$cluster[pop_info$siteID=="11_KLM"] <- "west"
pop_info$cluster[pop_info$siteID=="12_MSN"] <- "west"
pop_info$cluster[pop_info$siteID=="13_MSS"] <- "west"
pop_info$cluster[pop_info$siteID=="14_NGK"] <- "west"
pop_info$cluster[pop_info$siteID=="15_NGU"] <- "west"
pop_info$cluster[pop_info$siteID=="16_MNP"] <- "east"
pop_info$cluster[pop_info$siteID=="17_HND"] <- "east"
pop_info$cluster[pop_info$siteID=="18_AMR"] <- "east"
pop_info$cluster[pop_info$siteID=="19_CNP"] <- "east"
pop_info$cluster[pop_info$siteID=="20_KIB"] <- "east"
pop_info$cluster[pop_info$siteID=="21_TSW"] <- "east"
pop_info$cluster[pop_info$siteID=="22_KIN"] <- "east"
pop_info$cluster[pop_info$siteID=="23_SHI"] <- "east"
pop_info$cluster[pop_info$siteID=="24_SHT"] <- "east"
myData@other$cluster <- pop_info$cluster



#remove some PCA and GPS coordinate outliers and split TSW into two "populations"
##############################################
#build myData5
removeInd <- c("KIB_Kibwe013","IKR_IKR002","MRB_MMR231","GTR_GTR008","KAP_KAP066","KAP_KAP009","MSN_MSA001","IKR_IKR013","SHT_SHti082","SHI_SHpe034","SHI_SHpe035","SHI_SHma039","SHI_SHma044","SHI_SHma045","SHI_SHma046","SHI_SHma038","SHI_SHma040")
myData5<- myData[!row.names(myData@tab) %in% removeInd]

#build myData6
removeInd2 <-c("KIB_Kibwe001","CNP_CNP098","CNP_CNP177","CNP_CNP176","MRT_MMR119","MRT_MMR127","RUM_Ruma027","RUM_Ruma028","RUM_Ruma032","RUM_Ruma035","RUM_Ruma029")
myData6 <- myData5[!row.names(myData5@tab) %in% removeInd2]
tswa <- c("TSW_TSWng013","TSW_TSWng019",
          "TSW_TSWng047","TSW_TSWng048","TSW_TSWng053","TSW_TSWng056","TSW_TSWng065","TSW_TSWng066","TSW_TSWng074","TSW_TSWng076","TSW_TSWng077","TSW_TSWng095","TSW_TSWng096","TSW_TSWng097","TSW_TSWng098")
tswb <- row.names(myData6@tab)[!row.names(myData6@tab) %in% tswa & myData6@pop=="21_TSW"]
levels(myData6@pop) <- c(levels(myData5@pop),"21a_TSWa", "21b_TSWb")
myData6@pop[row.names(myData6@tab) %in% tswa] <- "21a_TSWa"
myData6@pop[row.names(myData6@tab) %in% tswb] <- "21b_TSWb"
myData6@pop <- droplevels(myData6@pop)
myData6@pop <- factor(myData6@pop,levels=c("01_KAP","02_OTT","03_RUM","04_GVR","05_MRB","06_MRT","07_NBS","08_FGT","09_IKR","10_GTR","11_KLM","12_MSN","13_MSS","14_NGK","15_NGU","16_MNP","17_HND","18_AMR","19_CNP","20_KIB","21a_TSWa","21b_TSWb","22_KIN","23_SHI","24_SHT"))

#build myData7 March 2, 2020
myData7 <- myData6
gvra <- c("GVR_MMR162","GVR_MMR164","GVR_MMR165","GVR_MMR166","GVR_MMR167","GVR_MMR191","GVR_MMR192","GVR_MMR194","GVR_MMR196","GVR_MMR197","GVR_MMR198")
gvrb <- row.names(myData7@tab)[!row.names(myData7@tab) %in% gvra & myData7@pop=="04_GVR"]
kiba <- c("KIB_Kibwe010","KIB_Kibwe011","KIB_Kibwe012","KIB_Kibwe002","KIB_Kibwe003","KIB_Kibwe004","KIB_Kibwe005","KIB_Kibwe006","KIB_Kibwe008")
kibb <- c("KIB_Kibwe025","KIB_Kibwe026","KIB_Kibwe024","KIB_Kibwe037","KIB_Kibwe038","KIB_Kibwe021","KIB_Kibwe022","KIB_Kibwe023","KIB_Kibwe016","KIB_Kibwe017","KIB_Kibwe018","KIB_Kibwe019","KIB_Kibwe020")
cnpa <- c("CNP_CNP039","CNP_CNP040","CNP_CNP041","CNP_CNP050","CNP_CNP051","CNP_CNP100","CNP_CNP101","CNP_CNP102","CNP_CNP103","CNP_CNP118","CNP_CNP119")
cnpb <- c("CNP_CNP001","CNP_CNP002","CNP_CNP003","CNP_CNP004","CNP_CNP005","CNP_CNP092","CNP_CNP093","CNP_CNP094","CNP_CNP095","CNP_CNP096")
kina <- c("KIN_KINny016","KIN_KINny017","KIN_KINny018","KIN_KINny019","KIN_KINny020","KIN_KINny021","KIN_KINny022","KIN_KINny024","KIN_KINny025","KIN_KINny027")
kinb <- c("KIN_KINny028","KIN_KINny029","KIN_KINny030","KIN_KINny031","KIN_KINny032","KIN_KINny033","KIN_KINny034","KIN_KINny035","KIN_KINny038","KIN_KINny042","KIN_KINny043","KIN_KINny044","KIN_KINny100")

levels(myData7@pop) <- c(levels(myData7@pop),"04a_GVRa", "04b_GVRb","19a_CNPa","19b_CNPb","20a_KIBa","20b_KIBb","22a_KINa","22b_KINb")
myData7@pop[row.names(myData7@tab) %in% gvra] <- "04a_GVRa"
myData7@pop[row.names(myData7@tab) %in% gvrb] <- "04b_GVRb"
myData7@pop[row.names(myData7@tab) %in% cnpa] <- "19a_CNPa"
myData7@pop[row.names(myData7@tab) %in% cnpb] <- "19b_CNPb"
myData7@pop[row.names(myData7@tab) %in% kiba] <- "20a_KIBa"
myData7@pop[row.names(myData7@tab) %in% kibb] <- "20b_KIBb"
myData7@pop[row.names(myData7@tab) %in% kina] <- "22a_KINa"
myData7@pop[row.names(myData7@tab) %in% kinb] <- "22b_KINb"
myData7@pop <- droplevels(myData7@pop)

myData8 <- myData7[myData7@pop %in% c("01_KAP","02_OTT","03_RUM","04a_GVRa","04b_GVRb","05_MRB","06_MRT","07_NBS","08_FGT","09_IKR","10_GTR","11_KLM","12_MSN","13_MSS","14_NGK","15_NGU","16_MNP","17_HND","18_AMR","19a_CNPa","19b_CNPb","20a_KIBa","20b_KIBb","21a_TSWa","21b_TSWb","22a_KINa","22b_KINb","23_SHI","24_SHT")]

myData8@pop <- factor(myData8@pop,levels=c("01_KAP","02_OTT","03_RUM","04a_GVRa","04b_GVRb","05_MRB","06_MRT","07_NBS","08_FGT","09_IKR","10_GTR","11_KLM","12_MSN","13_MSS","14_NGK","15_NGU","16_MNP","17_HND","18_AMR","19a_CNPa","19b_CNPb","20a_KIBa","20b_KIBb","21a_TSWa","21b_TSWb","22a_KINa","22b_KINb","23_SHI","24_SHT"))



#calculate GPS centroids
##############################################
siteID  <- myData8@pop

#define what you want to plot and create a dataframe:
corLat <- myData8@other$e_xy[,c("corLat")]
corLong <- myData8@other$e_xy[,c("corLong")]

df <- data.frame(siteID,corLong,corLat)
GPScentroids <- aggregate(cbind(corLong,corLat)~siteID,df,mean)
names(GPScentroids) <- c("siteID","centroidLong","centroidLat")
pop_info_centroids <- left_join(df,GPScentroids,by="siteID")

#place centroids into myData8 object and write new pop info to file
myData8@other$c_xy <- pop_info_centroids[,c("centroidLong","centroidLat")]



#assign closest pixel coordinates to ultimately save on compute time
##############################################
library(geostatsp)

#Load a stack of raster layers that have the pixels we will be using:
rast_stack <- stack("chelsa_merit_vars_kenya.tif")

#isolate one raster layer 
rast <- as(subset(rast_stack, 1),"RasterLayer")

#convert the raster layer to class "im" (image)
rast_im <- as.im.RasterLayer(rast)

#find the nearest pixels to exact coordinates
nearest_pixels_e <- nearest.raster.point(myData8@other$e_xy$corLong,myData8@other$e_xy$corLat,rast_im, indices=FALSE)
npe <- data.frame(nearest_pixels_e$x,nearest_pixels_e$y)
names(npe) <- c("pixelLong","pixelLat")

#find the nearest pixels to centroid coordinates
nearest_pixels_c <- nearest.raster.point(myData8@other$c_xy$centroidLong,myData8@other$c_xy$centroidLat,rast_im, indices=FALSE)
npc <- data.frame(nearest_pixels_c$x,nearest_pixels_c$y)
names(npc) <- c("pixelcentroidLong","pixelcentroidLat")

pop_info_pixels <- cbind(rownames(myData8@tab),myData8@other$region,myData8@other$cluster,npe,npc)
names(pop_info_pixels) <- c("indivID","region","cluster","pixelLong","pixelLat","pixelcentroidLong","pixelcentroidLat")

#place pixel coordinates into myData8 object
myData8@other$pix_e_xy <- npe[,c("pixelLong","pixelLat")]
myData8@other$pix_c_xy <- npc[,c("pixelcentroidLong","pixelcentroidLat")]


#write the final myData8 genepop and info to file so it can be loaded directly
##############################################
pop_info_centroids_pixels <- cbind(pop_info_pixels,pop_info_centroids)
write.csv(pop_info_centroids_pixels[order(pop_info_centroids_pixels$siteID),c("indivID","siteID","pixelLat","pixelLong","corLat","corLong","centroidLat","centroidLong","region","cluster")],file="Gpd_KenTza_11loci_659indv_NPSmyData8_info.csv")
library(graph4lg)
genind_to_genepop(myData8,output="Gpd_KenTza_11loci_659indv_NPSmyData8.txt")
#edit to .gen to match format that I like with one line for title and locus names and remove siteID from indivID

#reload the final data as myData9 here to be sure we are working with exactly the right siteID and indvID
##############################################
myData9 <- read.genepop("Gpd_KenTza_11loci_659indv_NPSmyData8.gen",ncode=3)
myData9_info <- read.table("Gpd_KenTza_11loci_659indv_NPSmyData8_info.csv",header=T,sep=",")

#add the name of the populations to the geneind object
myData9@pop <- myData9_info$siteID
myData9@other$e_xy <- myData9_info[,c("corLong","corLat")]
myData9@other$pixel_e_xy <- myData9_info[,c("pixelLong","pixelLat")]
myData9@other$c_xy <- myData9_info[,c("centroidLong","centroidLat")]
myData9@other$region <- myData9_info$region
myData9@other$cluster <- myData9_info$cluster


#############################################################################
#Complete DAPC with adegenet following tutorial at http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
#############################################################################
#DAPC 
##############################################
#interpret number of clusters
#bic <- find.clusters(myData,max.n.clust=23,n.pca = 500) #keep all PCs, just put 2 clusters so you can see the plot
grp2 <- find.clusters(myData9,n.pca=4,n.clust=2) 
grp3 <- find.clusters(myData9,n.pca=4,n.clust=3) 
grp4 <- find.clusters(myData9,n.pca=4,n.clust=4) 
grp5 <- find.clusters(myData9,n.pca=4,n.clust=5) 
grp6 <- find.clusters(myData9,n.pca=4,n.clust=6) 
grp7 <- find.clusters(myData9,n.pca=4,n.clust=7) 
grp8 <- find.clusters(myData9,n.pca=4,n.clust=8) 

#Group by DAPC results at different K values
grp2 <-grp2$grp
dapc2 <- dapc(myData9, grp2,n.pca=4,n.da=2) 
grp3 <-grp3$grp
dapc3 <- dapc(myData9, grp3,n.pca=4,n.da=2) 
grp4 <-grp4$grp
dapc4 <- dapc(myData9, grp4,n.pca=4,n.da=2) 
grp5 <-grp5$grp
dapc5 <- dapc(myData9, grp5,n.pca=4,n.da=2) 
grp6 <-grp6$grp
dapc6 <- dapc(myData9, grp6,n.pca=4,n.da=2) 
grp7 <-grp7$grp
dapc7 <- dapc(myData9, grp7,n.pca=4,n.da=2) 
grp8 <-grp8$grp
dapc8 <- dapc(myData9, grp8,n.pca=4,n.da=2) 

par(mfrow = c(1, 1))
table.value(table(pop(myData9), grp2), col.lab=paste("Group", 1:2), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)
table.value(table(pop(myData9), grp3), col.lab=paste("Group", 1:3), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)
table.value(table(pop(myData9), grp4), col.lab=paste("Group", 1:4), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)
table.value(table(pop(myData9), grp5), col.lab=paste("Group", 1:5), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)
table.value(table(pop(myData9), grp6), col.lab=paste("Group", 1:6), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)
table.value(table(pop(myData9), grp7), col.lab=paste("Group", 1:7), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)
table.value(table(pop(myData9), grp8), col.lab=paste("Group", 1:8), row.lab=paste(sort(unique(myData9@pop))), clabel.row = .75, clegend = 0)

par(mfrow=c(3,2),mar=c(0,5,5,1)+.1)
col.green <- "#0C660C"
col.pink <-"#980043"
col.orange <- "#FF6600"
col.blue <- "#0099E6"

#K2 
scatter(dapc2, scree.da=FALSE,scree.pca = FALSE, posi.pca = "bottomright", ratio.pca=.25,col=c(col.blue,col.green))
compoplot(dapc2, legend=FALSE,posi="topright",txt.leg=paste("Cluster", 1:2), lab=myData9@pop,ncol=2, xlab="individuals", col=c(col.blue,col.green),show.lab=T, subset=row.names(myData9@tab)[order(myData9@pop)])

#K3
scatter(dapc3, scree.da=FALSE,scree.pca = FALSE, posi.pca = "bottomright", ratio.pca=.25, col=c(col.blue,col.green,col.orange))
compoplot(dapc3, legend=FALSE,posi="topright",txt.leg=paste("Cluster", 1:3), lab=myData9@pop,ncol=2, xlab="individuals", col=c(col.blue,col.green,col.orange),show.lab=T, subset=row.names(myData9@tab)[order(myData9@pop)])

#K4
scatter(dapc4, scree.da=FALSE,scree.pca = FALSE,posi.pca = "bottomright", ratio.pca=.25, col=c(col.blue,col.green,col.orange,col.pink))
compoplot(dapc4, legend=FALSE,posi="topright",txt.leg=paste("Cluster", 1:4), lab=myData9@pop,ncol=2, xlab="individuals", col=c(col.blue,col.green,col.orange,col.pink),show.lab=T, subset=row.names(myData9@tab)[order(myData9@pop)])

#K5
scatter(dapc5, scree.da=FALSE,scree.pca = FALSE, posi.pca = "bottomright", ratio.pca=.25, col=c("#33004C","#0099E6","#FF6600","#99194C","#0C660C"))
compoplot(dapc5, legend=FALSE,posi="topright",txt.leg=paste("Cluster", 1:5), lab=myData9@pop,ncol=2, xlab="individuals", col=c("#33004C","#0099E6","#FF6600","#99194C","#0C660C"),show.lab=T, subset=row.names(myData9@tab)[order(myData9@pop)])

#K6
scatter(dapc6, scree.da=FALSE,scree.pca = FALSE, posi.pca = "bottomright", ratio.pca=.25, col=c("#33004C","#FE9800","#FF6600","#99194C","#0099E6", "#0C660C"))
compoplot(dapc6, legend=FALSE,posi="topright",txt.leg=paste("Cluster", 1:6), lab=myData9@pop,ncol=2, xlab="individuals", col=c("#33004C","#FE9800","#FF6600","#99194C","#0099E6", "#0C660C"),show.lab=T, subset=row.names(myData9@tab)[order(myData9@pop)])

#K7
scatter(dapc7, scree.da=FALSE,scree.pca = FALSE, posi.pca = "bottomright", ratio.pca=.25, 
        col=c("#FE9800","#33004C","#0099E6","#856f2c","#0C660C", "#99194C","#FF6600"))
compoplot(dapc7,legend=FALSE, posi="topright",txt.leg=paste("Cluster", 1:7), lab=myData9@pop,ncol=2, xlab="individuals", col=c("#FE9800","#33004C","#0099E6","#856f2c","#0C660C", "#99194C","#FF6600"),show.lab=T, subset=row.names(myData9@tab)[order(myData9@pop)])


#create pairwise distance table by individual
##############################################
library(philentropy)

################################
#calculate PCA distance with 10 PCs
################################
#There are missing data. They will all be replaced by scaleGen:
X <- scaleGen(myData9, NA.method="mean")
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=30)
#find percent variance explained for each of the three PCs
pc1 <- round(pca1$eig[1]/sum(pca1$eig)*100,digits=2)
pc2 <- round(pca1$eig[2]/sum(pca1$eig)*100,digits=2)
pc3 <- round(pca1$eig[3]/sum(pca1$eig)*100,digits=2)
#assign PC axis to plot
PC1 <- pca1$li$Axis1
PC2 <- pca1$li$Axis2
PC3 <- pca1$li$Axis3
#view eigenvalues
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

##################################
#visualize by region
##################################
cluster <- myData9@other$region
df <- data.frame(cluster,PC1,PC2,PC3)
centroids <- aggregate(cbind(PC1,PC2,PC3)~cluster,df,mean)
gg <- merge(df,aggregate(cbind(mean.x=PC1,mean.y=PC2,mean.z=PC3)~cluster,df,mean),by="cluster")
psnp <- ggplot(data = gg, aes(PC1,PC2,color=factor(cluster)),)+
  xlab(paste("PC 1 (",pc1,"%)"))+
  ylab(paste("PC 2 (",pc2,"%)"))+
  geom_point(size=2)+
  #stat_ellipse(type="norm", level = 0.9, linetype=1)+
  geom_point(aes(x=mean.x,y=mean.y),size=0)+
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2))+
  scale_colour_manual(values = c(rainbow((length(unique(cluster)))+1)))+
  ggtitle("PCA by geographic region")+
  theme(plot.title = element_text(hjust = 0.5))+
  #scale_x_reverse()+
  theme(panel.grid.major = element_line(colour = "#856f2c"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.key = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
psnp
ggsave(plot=psnp,height=5,width=6,dpi=200, filename="Gpd_PCA_byRegion.pdf", useDingbats=FALSE)
plot_ly(x=PC1, y=PC2, z=PC3, type="scatter3d", mode="markers", color=factor(cluster),colors=rainbow((length(unique(cluster)))+1))
s.label(pca1$li)

##################################
#visualize by sampling site
##################################
cluster <- myData9@pop
df <- data.frame(cluster,PC1,PC2,PC3)
centroids <- aggregate(cbind(PC1,PC2,PC3)~cluster,df,mean)
gg <- merge(df,aggregate(cbind(mean.x=PC1,mean.y=PC2,mean.z=PC3)~cluster,df,mean),by="cluster")
psnp <- ggplot(data = gg, aes(PC1,PC2,color=factor(cluster)),)+
  xlab(paste("PC 1 (",pc1,"%)"))+
  ylab(paste("PC 2 (",pc2,"%)"))+
  geom_point(size=2)+
  #stat_ellipse(type="norm", level = 0.9, linetype=1)+
  geom_point(aes(x=mean.x,y=mean.y),size=0)+
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2))+
  scale_colour_manual(values = c(rainbow((length(unique(cluster)))+1)))+
  ggtitle("PCA by sampling site")+
  theme(plot.title = element_text(hjust = 0.5))+
  #scale_x_reverse()+
  theme(panel.grid.major = element_line(colour = "#856f2c"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA), legend.key = element_blank(), legend.title = element_blank(), axis.line = element_line(colour = "black"), aspect.ratio = 1)
psnp
ggsave(plot=psnp,height=5,width=6,dpi=200, filename="Gpd_PCA_bySite.pdf", useDingbats=FALSE)
plot_ly(x=PC1, y=PC2, z=PC3, type="scatter3d", mode="markers", color=factor(cluster),colors=rainbow((length(unique(cluster)))+1))
s.label(pca1$li)

##################################
#calculate euclidean distance with first 10 PCs
##################################
pca_matrix <- distance(pca1$li[,c(1:10)], method = "euclidean")
pca_upper <- pca_matrix
colnames(pca_upper) <- rownames(pca1$li)
rownames(pca_upper) <- rownames(pca1$li)
pca_upper[lower.tri(pca_upper,diag=TRUE)] <- NA

#convert into a two column table instead of a square matrix
ind_pairwise <- na.omit(melt(pca_upper))
names(ind_pairwise) <- c("Ind1","Ind2","dist10PCs")


#calculate geodistance by individual
##############################################
#join ind_pairwise with geographic data
ind_data <- cbind(rownames(myData9@tab)[order(myData9@pop)],myData9@pop[order(myData9@pop)],myData9@other$e_xy[order(myData9@pop),1:2],myData9@other$pixel_e_xy[order(myData9@pop),1:2],myData9@other$region[order(myData9@pop)],myData9@other$cluster[order(myData9@pop)])
names(ind_data) <- c("indivID","siteID","corLong","corLat","pixelLong","pixelLat","region","cluster")
ind_uniq <- distinct(ind_data)

Ind1_join <- left_join(ind_pairwise, ind_data, by = c("Ind1" = "indivID"))
Ind2_join <- left_join(ind_pairwise, ind_data, by = c("Ind2" = "indivID"))
names(Ind1_join) <- c(names(Ind1_join)[1:3],paste0("Ind1_",names(Ind1_join)[-1:-3]))
names(Ind2_join) <- c(names(Ind2_join)[1:3],paste0("Ind2_",names(Ind2_join)[-1:-3]))
Ind1_Ind2_join <- left_join(Ind1_join,Ind2_join)

indiv_pairwise <- distinct(Ind1_Ind2_join)

begin.table <- indiv_pairwise[,c("Ind1_pixelLong","Ind1_pixelLat")]
begin.coord <- begin.table
begin.Sp.points <- SpatialPoints(begin.coord)

end.table <- indiv_pairwise[,c("Ind2_pixelLong","Ind2_pixelLat")]
end.coord <- end.table
end.Sp.points <- SpatialPoints(end.coord)

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))
spatial.p <- as(p, "SpatialLines")
proj4string(spatial.p) <- crs.geo  # define projection system of our data
geodist <- pointDistance(begin.Sp.points, end.Sp.points, lonlat=TRUE, allpairs=FALSE)

#add m and km geodist to dataframe
indiv_pairwise$m <- geodist
indiv_pairwise$km <- indiv_pairwise$m/1000

#save the indiv_pairwise as a table
write.table(indiv_pairwise,file="Gpd_KenTza_11loci_659indv_indiv_pairwise.txt")


#quantify within-site distances
##############################################
within_sites  <- na.omit(indiv_pairwise[indiv_pairwise$Ind1_siteID==indiv_pairwise$Ind2_siteID,])

#make a histogram:
hist(within_sites$km[within_sites$Ind1!=within_sites$Ind2],breaks=max(within_sites$km),main="Within sampling site distances \n(based on pixel coordinates)")

#what is the percent of pairs > 3 km apart per site?
i=0
sdf <- data.frame("siteID"=sort(unique(within_sites$Ind1_siteID)),"percent3"=0)
for (i in 1:length(sdf$siteID))
{
  sdf$percent3[i] <- (length(within_sites$km[within_sites$Ind1_siteID==sdf$siteID[i] & within_sites$km>3 & within_sites$Ind1!=within_sites$Ind2])/length(within_sites$km[within_sites$Ind1_siteID==sdf$siteID[i] & within_sites$Ind1!=within_sites$Ind2]))
}
sdt <- data.frame("siteID"="total","percent3"=length(within_sites$km[within_sites$km>3 & within_sites$Ind1!=within_sites$Ind2])/length(within_sites$km[within_sites$Ind1!=within_sites$Ind2]))
sdf <- rbind(sdf,sdt)

#what is the percent of pairs > 2 km apart per site?
i=0
sdf2 <- data.frame("siteID"=sort(unique(within_sites$Ind1_siteID)),"percent2"=0)
for (i in 1:length(sdf2$siteID))
{
  sdf2$percent2[i] <- (length(within_sites$km[within_sites$Ind1_siteID==sdf2$siteID[i] & within_sites$km>2 & within_sites$Ind1!=within_sites$Ind2])/length(within_sites$km[within_sites$Ind1_siteID==sdf2$siteID[i] & within_sites$Ind1!=within_sites$Ind2]))
}
sdt2 <- data.frame("siteID"="total","percent2"=length(within_sites$km[within_sites$km>2 & within_sites$Ind1!=within_sites$Ind2])/length(within_sites$km[within_sites$Ind1!=within_sites$Ind2]))
sdf2 <- rbind(sdf2,sdt2)
sdf <- cbind(sdf,sdf2)

#barplot
barplot(sdf$percent3,names.arg=sdf$siteID,las=3,horiz=F,cex.names=.6,main="Proportion of within-site distances > 3 km")
barplot(sdf$percent2,names.arg=sdf$siteID,las=3,horiz=F,cex.names=.6,main="Proportion of within-site distances > 2 km")

sdf


#plot trap pixel coordinates for areas of concern, include=FALSE
##############################################
ind_data <- cbind(rownames(myData9@tab)[order(myData9@pop)],myData9@pop[order(myData9@pop)],myData9@other$e_xy[order(myData9@pop),1:2],myData9@other$pixel_e_xy[order(myData9@pop),1:2],myData9@other$region[order(myData9@pop)],myData9@other$cluster[order(myData9@pop)])
names(ind_data) <- c("indivID","siteID","corLong","corLat","pixelLong","pixelLat","region","cluster")
ind_uniq <- distinct(ind_data)

#CNP for example
ex_Points <- ind_uniq[ind_uniq$siteID %in% c("19a_CNPa"),]
ex_spPoints <-   SpatialPoints(ex_Points[,c("pixelLong","pixelLat")])
proj4string(ex_spPoints) <- crs.geo
map(,ylim=c((min(ex_Points$pixelLat)-.1),(max(ex_Points$pixelLat)+.1)),xlim=c((min(ex_Points$pixelLong)-.1),(max(ex_Points$pixelLong)+.1)))
points(ex_spPoints,pch=3)
text(ex_spPoints,labels=ex_Points$indivID,cex=.5,adj=c(-.5,0),col="blue") 

ex_Points <- ind_uniq[ind_uniq$siteID %in% c("19b_CNPb"),]
ex_spPoints <-   SpatialPoints(ex_Points[,c("pixelLong","pixelLat")])
proj4string(ex_spPoints) <- crs.geo
map(,ylim=c((min(ex_Points$pixelLat)-.1),(max(ex_Points$pixelLat)+.1)),xlim=c((min(ex_Points$pixelLong)-.1),(max(ex_Points$pixelLong)+.1)))
points(ex_spPoints,pch=3)
text(ex_spPoints,labels=ex_Points$indivID,cex=.5,adj=c(-.5,0),col="blue") 

ex_Points <- ind_uniq[ind_uniq$siteID %in% c("19a_CNPa"),]
ex_spPoints <-   SpatialPoints(ex_Points[,c("corLong","corLat")])
proj4string(ex_spPoints) <- crs.geo
map(,ylim=c((min(ex_Points$pixelLat)-.1),(max(ex_Points$pixelLat)+.1)),xlim=c((min(ex_Points$pixelLong)-.1),(max(ex_Points$pixelLong)+.1)))
points(ex_spPoints,pch=3)
text(ex_spPoints,labels=ex_Points$indivID,cex=.5,adj=c(-.5,0),col="blue") 

ex_Points <- ind_uniq[ind_uniq$siteID %in% c("19b_CNPb"),]
ex_spPoints <-   SpatialPoints(ex_Points[,c("corLong","corLat")])
proj4string(ex_spPoints) <- crs.geo
map(,ylim=c((min(ex_Points$pixelLat)-.1),(max(ex_Points$pixelLat)+.1)),xlim=c((min(ex_Points$pixelLong)-.1),(max(ex_Points$pixelLong)+.1)))
points(ex_spPoints,pch=3)
text(ex_spPoints,labels=ex_Points$indivID,cex=.5,adj=c(-.5,0),col="blue") 



#map centroids and within-site straight lines between pixel localities
##############################################
coords <- data.frame(cbind(myData9@pop[order(myData9@pop)],myData9@other$c_xy[order(myData9@pop),]))
unique_coords <- distinct(coords)
names(unique_coords) <- c("siteID","centroidLong","centroidLat")
spPoints <- SpatialPoints(unique_coords[,c("centroidLong","centroidLat")])

#make spatial lines to plot:
begin.table <- within_sites[,c("Ind1_pixelLong","Ind1_pixelLat")]
begin.coord <- begin.table
begin.Sp.points <- SpatialPoints(begin.coord)
end.table <- within_sites[,c("Ind2_pixelLong","Ind2_pixelLat")]
end.coord <- end.table
end.Sp.points <- SpatialPoints(end.coord)
p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))
spatial.p <- as(p, "SpatialLines")
proj4string(spatial.p) <- crs.geo  # define projection system of our data

#plot
map(,ylim=c(-5,6),xlim=c(33.5,43))
points(spPoints,pch=19,col=c(rainbow(length(unique_coords$siteID)+1)))
text(spPoints,labels=unique_coords$siteID,cex=.5,adj=c(-.2,0))
plot(spatial.p,add=T,main="Within sites straight line distances")

#save to pdf:
pdf(file="Gpd_localities_within_site_straight_lines_check.pdf",paper="letter")
map(,ylim=c(-5,6),xlim=c(33.5,43))
points(spPoints,pch=19,col=c(rainbow(length(unique_coords$siteID)+1)))
text(spPoints,labels=unique_coords$siteID,cex=.5,adj=c(-.2,0))
plot(spatial.p,add=T,main="Within sites straight line distances")
dev.off()


###################################
#IBD by individual
###################################
#isolation by distance by individual
##############################################
G.table <- indiv_pairwise

#cannot have distance of 0 or 1 or will create difficulty in IBD
G.table$m[G.table$Ind1_siteID==G.table$Ind2_siteID] <- 500
G.table$km[G.table$Ind1_siteID==G.table$Ind2_siteID] <- .5

Dgen_df <- as.data.frame(G.table[,c("Ind1","Ind2","dist10PCs")])
Dgen_wide <- reshape(Dgen_df,v.names="dist10PCs",idvar="Ind1",timevar="Ind2",direction="wide")
rownames(Dgen_wide) <- Dgen_wide[,1]
Dgen <- as.dist(t(Dgen_wide[,-1]))

Dgeo_df <- as.data.frame(G.table[,c("Ind1","Ind2","km")])
Dgeo_wide <- reshape(Dgeo_df,v.names="km",idvar="Ind1",timevar="Ind2",direction="wide")
rownames(Dgeo_wide) <- Dgeo_wide[,1]
Dgeo <- as.dist(t(Dgeo_wide[,-1]))

all.ibd <- mantel.randtest(Dgen,Dgeo)
all.ibd
all.lm <- lm(Dgen_df$dist10PCs~Dgeo_df$km)
all.lm

par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(all.ibd)
plot(Dgeo_df$km,Dgen_df$dist10PCs,xlab="km",ylab="PCA distance (10 PCs)",main="All pairwise")
abline(all.lm,col="red")


#isolation by distance for east only by individual
##############################################
east.G.table <- indiv_pairwise[indiv_pairwise$Ind1_cluster=="east" & indiv_pairwise$Ind2_cluster=="east",]

#cannot have distance of 0 or 1 or will create difficulty in IBD
east.G.table$m[east.G.table$Ind1_siteID==east.G.table$Ind2_siteID] <- 500
east.G.table$km[east.G.table$Ind1_siteID==east.G.table$Ind2_siteID] <- .5

Dgen_df <- as.data.frame(east.G.table[,c("Ind1","Ind2","dist10PCs")])
Dgen_wide <- reshape(Dgen_df,v.names="dist10PCs",idvar="Ind1",timevar="Ind2",direction="wide")
rownames(Dgen_wide) <- Dgen_wide[,1]
Dgen <- as.dist(t(Dgen_wide[,-1]))

Dgeo_df <- as.data.frame(east.G.table[,c("Ind1","Ind2","km")])
Dgeo_wide <- reshape(Dgeo_df,v.names="km",idvar="Ind1",timevar="Ind2",direction="wide")
rownames(Dgeo_wide) <- Dgeo_wide[,1]
Dgeo <- as.dist(t(Dgeo_wide[,-1]))

east.ibd <- mantel.randtest(Dgen,Dgeo)
east.ibd
east.lm <- lm(Dgen_df$dist10PCs~Dgeo_df$km)
east.lm

par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(east.ibd)
plot(Dgeo_df$km,Dgen_df$dist10PCs,xlab="km",ylab="PCA distance (10 PCs)",main="East only")
abline(east.lm,col="red")


#isolation by distance for west only by individual
##############################################
west.G.table <- indiv_pairwise[indiv_pairwise$Ind1_cluster=="west" & indiv_pairwise$Ind2_cluster=="west",]

#cannot have distance of 0 or 1 or will create difficulty in IBD
west.G.table$m[west.G.table$Ind1_siteID==west.G.table$Ind2_siteID] <- 500
west.G.table$km[west.G.table$Ind1_siteID==west.G.table$Ind2_siteID] <- .5

Dgen_df <- as.data.frame(west.G.table[,c("Ind1","Ind2","dist10PCs")])
Dgen_wide <- reshape(Dgen_df,v.names="dist10PCs",idvar="Ind1",timevar="Ind2",direction="wide")
rownames(Dgen_wide) <- Dgen_wide[,1]
Dgen <- as.dist(t(Dgen_wide[,-1]))

Dgeo_df <- as.data.frame(west.G.table[,c("Ind1","Ind2","km")])
Dgeo_wide <- reshape(Dgeo_df,v.names="km",idvar="Ind1",timevar="Ind2",direction="wide")
rownames(Dgeo_wide) <- Dgeo_wide[,1]
Dgeo <- as.dist(t(Dgeo_wide[,-1]))

west.ibd <- mantel.randtest(Dgen,Dgeo)
west.ibd
west.lm <- lm(Dgen_df$dist10PCs~Dgeo_df$km)
west.lm

par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(west.ibd)
plot(Dgeo_df$km,Dgen_df$dist10PCs,xlab="km",ylab="PCA distance (10 PCs)",main="West only")
abline(west.lm,col="red")


###################################
#IBD by population
###################################
#reload the final data here to be sure we are working with exactly the right siteID and indvID
##############################################
myData9 <- read.genepop("Gpd_KenTza_11loci_659indv_NPSmyData8.gen",ncode=3)
myData9_info <- read.table("Gpd_KenTza_11loci_659indv_NPSmyData8_info.csv",header=T,sep=",")

#add the name of the populations to the geneind object
myData9@pop <- myData9_info$siteID
myData9@other$e_xy <- myData9_info[,c("corLong","corLat")]
myData9@other$pixel_e_xy <- myData9_info[,c("pixelLong","pixelLat")]
myData9@other$c_xy <- myData9_info[,c("centroidLong","centroidLat")]
myData9@other$region <- myData9_info$region
myData9@other$cluster <- myData9_info$cluster



#build a pop_pairwise table
##############################################
################################
#calculate pairwise genetic distance by sampling site 
################################
pop_fst_matrix <- pairwise.fst(myData9,res.type="matrix")
pop_fst_upper <- pop_fst_matrix
pop_fst_upper[lower.tri(pop_fst_upper,diag=TRUE)] <- NA
pop_Linearfst_matrix <- (pop_fst_matrix)/(1-(pop_fst_matrix))
#convert into a two column table instead of a square matrix
pop_fst_table <- na.omit(melt(pop_fst_upper))
names(pop_fst_table) <- c("Pop1","Pop2","Fst")
#Add in linearized Fst based on Rousset's equation
pop_fst_table$LinearFst <- (pop_fst_table$Fst)/(1-(pop_fst_table$Fst))

################################
#calculate Reynold's distance
################################
myPop <- genind2genpop(myData9)
pop_Rdist_matrix <- as.matrix(dist.genpop(myPop, method=3,upper=T,diag=T))
pop_Rdist_upper <- pop_Rdist_matrix
pop_Rdist_upper[lower.tri(pop_Rdist_upper,diag=TRUE)] <- NA
#convert into a two column table instead of a square matrix
pop_Rdist_table <- na.omit(melt(pop_Rdist_upper))
names(pop_Rdist_table) <- c("Pop1","Pop2","Reynolds")

################################
#calculate Edward's distance
################################
pop_Edist_matrix <- as.matrix(dist.genpop(myPop, method=2,upper=T,diag=T))
pop_Edist_upper <- pop_Edist_matrix
pop_Edist_upper[lower.tri(pop_Edist_upper,diag=TRUE)] <- NA
#convert into a two column table instead of a square matrix
pop_Edist_table <- na.omit(melt(pop_Edist_upper))
names(pop_Edist_table) <- c("Pop1","Pop2","Edwards")

################################
#Create pairwise genetic distance table
################################
pop_pairwise <-  pop_fst_table
pop_pairwise$Reynolds <- pop_Rdist_table$Reynolds
pop_pairwise$Edwards <- pop_Edist_table$Edwards


#calculate geodistance by pop
##############################################
#join pop_pairwise with geographic data
pop_data <- cbind(myData9@pop[order(myData9@pop)],myData9@other$c_xy[order(myData9@pop),1:2],myData9@other$region[order(myData9@pop)],myData9@other$cluster[order(myData9@pop)])
names(pop_data) <- c("siteID","centroidLong","centroidLat","region","cluster")
pop_uniq <- distinct(pop_data)

Pop1_join <- left_join(pop_pairwise, pop_uniq, by = c("Pop1" = "siteID"))
Pop2_join <- left_join(pop_pairwise, pop_uniq, by = c("Pop2" = "siteID"))
names(Pop1_join) <- c(names(Pop1_join)[1:6],paste0("Pop1_",names(Pop1_join)[-1:-6]))
names(Pop2_join) <- c(names(Pop2_join)[1:6],paste0("Pop2_",names(Pop2_join)[-1:-6]))
Pop1_Pop2_join <- left_join(Pop1_join,Pop2_join)

pop_pairwise <- distinct(Pop1_Pop2_join)

begin.table <- pop_pairwise[,c("Pop1_centroidLong","Pop1_centroidLat")]
begin.coord <- begin.table
begin.Sp.points <- SpatialPoints(begin.coord)

end.table <- pop_pairwise[,c("Pop2_centroidLong","Pop2_centroidLat")]
end.coord <- end.table
end.Sp.points <- SpatialPoints(end.coord)

p <- psp(begin.table[,1], begin.table[,2], end.table[,1], end.table[,2], owin(range(c(begin.table[,1], end.table[,1])), range(c(begin.table[,2], end.table[,2]))))
spatial.p <- as(p, "SpatialLines")
proj4string(spatial.p) <- crs.geo  # define projection system of our data
geodist <- pointDistance(begin.Sp.points, end.Sp.points, lonlat=TRUE, allpairs=FALSE)

#add m and km geodist to dataframe
pop_pairwise$m <- geodist
pop_pairwise$km <- pop_pairwise$m/1000

#write pop_pairwise to file as a table:
write.table(pop_pairwise,file="Gpd_KenTza_11loci_659indv_pop_pairwise.txt")

pop_pairwise$kmL <- log(pop_pairwise$m/1000)
pop_pairwise[pop_pairwise$km<5,]



#isolation by distance at population level
##############################################
G.table <- pop_pairwise
min(G.table$km)

Dgen_df <- as.data.frame(G.table[,c("Pop1","Pop2","LinearFst")])
Dgen_wide <- reshape(Dgen_df,v.names="LinearFst",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(Dgen_wide) <- Dgen_wide[,1]
Dgen <- as.dist(t(Dgen_wide[,-1]))

Dgeo_df <- as.data.frame(G.table[,c("Pop1","Pop2","km")])
Dgeo_wide <- reshape(Dgeo_df,v.names="km",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(Dgeo_wide) <- Dgeo_wide[,1]
Dgeo <- as.dist(t(Dgeo_wide[,-1]))

#log of geographic distance
DgeoL_df <- as.data.frame(G.table[,c("Pop1","Pop2","kmL")])
DgeoL_wide <- reshape(DgeoL_df,v.names="kmL",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(DgeoL_wide) <- DgeoL_wide[,1]
DgeoL <- as.dist(t(DgeoL_wide[,-1]))

ibd <- mantel.randtest(Dgen,Dgeo)
ibd
par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(ibd)
plot(Dgeo_df$km,Dgen_df$LinearFst,xlab="km",ylab="Linearized FST",main="All pairwise")
abline(lm(Dgen_df$LinearFst~Dgeo_df$km),col="red")

ibdL <- mantel.randtest(Dgen,DgeoL)
ibdL
par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(ibdL)
plot(DgeoL_df$kmL,Dgen_df$LinearFst,xlab="ln(km)",ylab="Linearized FST",main="All pairwise")
abline(lm(Dgen_df$LinearFst~DgeoL_df$kmL),col="red")

Dgen_melt <- na.omit(melt(Dgen_wide[,-1]))
Dgeo_melt <- na.omit(melt(Dgeo_wide[,-1]))
DgeoL_melt <- na.omit(melt(DgeoL_wide[,-1]))

lm_ibd <- lm(Dgen_melt$value~Dgeo_melt$value)
summary(lm_ibd)
lm_ibdL <- lm(Dgen_melt$value~DgeoL_melt$value)
summary(lm_ibdL)


#isolation by distance for east only
##############################################
east.p.G.table <- pop_pairwise[pop_pairwise$Pop1_cluster=="east" & pop_pairwise$Pop2_cluster=="east",]
min(east.p.G.table$km)

east.p.Dgen_df <- as.data.frame(east.p.G.table[,c("Pop1","Pop2","LinearFst")])
east.p.Dgen_wide <- reshape(east.p.Dgen_df,v.names="LinearFst",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(east.p.Dgen_wide) <- east.p.Dgen_wide[,1]
east.p.Dgen <- as.dist(t(east.p.Dgen_wide[,-1]))

east.p.Dgeo_df <- as.data.frame(east.p.G.table[,c("Pop1","Pop2","km")])
east.p.Dgeo_wide <- reshape(east.p.Dgeo_df,v.names="km",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(east.p.Dgeo_wide) <- east.p.Dgeo_wide[,1]
east.p.Dgeo <- as.dist(t(east.p.Dgeo_wide[,-1]))

#log of geographic distance
east.p.DgeoL_df <- as.data.frame(east.p.G.table[,c("Pop1","Pop2","kmL")])
east.p.DgeoL_wide <- reshape(east.p.DgeoL_df,v.names="kmL",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(east.p.DgeoL_wide) <- east.p.DgeoL_wide[,1]
east.p.DgeoL <- as.dist(t(east.p.DgeoL_wide[,-1]))

east.p.ibd <- mantel.randtest(east.p.Dgen,east.p.Dgeo)
east.p.ibd
par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(east.p.ibd)
plot(east.p.Dgeo_df$km,east.p.Dgen_df$LinearFst,xlab="km",ylab="Linearized FST",main="East only")
abline(lm(east.p.Dgen_df$LinearFst~east.p.Dgeo_df$km),col="red")

east.p.ibdL <- mantel.randtest(east.p.Dgen,east.p.DgeoL)
east.p.ibdL
par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(east.p.ibdL)
plot(east.p.DgeoL_df$kmL,east.p.Dgen_df$LinearFst,xlab="ln(km)",ylab="Linearized FST",main="East only")
abline(lm(east.p.Dgen_df$LinearFst~east.p.DgeoL_df$kmL),col="red")

east.p.Dgen_melt <- na.omit(melt(east.p.Dgen_wide[,-1]))
east.p.Dgeo_melt <- na.omit(melt(east.p.Dgeo_wide[,-1]))
east.p.DgeoL_melt <- na.omit(melt(east.p.DgeoL_wide[,-1]))

east.p.lm_ibd <- lm(east.p.Dgen_melt$value~east.p.Dgeo_melt$value)
summary(east.p.lm_ibd)
east.p.lm_ibdL <- lm(east.p.Dgen_melt$value~east.p.DgeoL_melt$value)
summary(east.p.lm_ibdL)


#isolation by distance for west only
##############################################
west.p.G.table <- pop_pairwise[pop_pairwise$Pop1_cluster=="west" & pop_pairwise$Pop2_cluster=="west",]
min(west.p.G.table$km)

west.p.Dgen_df <- as.data.frame(west.p.G.table[,c("Pop1","Pop2","LinearFst")])
west.p.Dgen_wide <- reshape(west.p.Dgen_df,v.names="LinearFst",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(west.p.Dgen_wide) <- west.p.Dgen_wide[,1]
west.p.Dgen <- as.dist(t(west.p.Dgen_wide[,-1]))

west.p.Dgeo_df <- as.data.frame(west.p.G.table[,c("Pop1","Pop2","km")])
west.p.Dgeo_wide <- reshape(west.p.Dgeo_df,v.names="km",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(west.p.Dgeo_wide) <- west.p.Dgeo_wide[,1]
west.p.Dgeo <- as.dist(t(west.p.Dgeo_wide[,-1]))

#log of geographic distance
west.p.DgeoL_df <- as.data.frame(west.p.G.table[,c("Pop1","Pop2","kmL")])
west.p.DgeoL_wide <- reshape(west.p.DgeoL_df,v.names="kmL",idvar="Pop1",timevar="Pop2",direction="wide")
rownames(west.p.DgeoL_wide) <- west.p.DgeoL_wide[,1]
west.p.DgeoL <- as.dist(t(west.p.DgeoL_wide[,-1]))

west.p.ibd <- mantel.randtest(west.p.Dgen,west.p.Dgeo)
west.p.ibd
par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(west.p.ibd)
plot(west.p.Dgeo_df$km,west.p.Dgen_df$LinearFst,xlab="km",ylab="Linearized FST",main="West only")
abline(lm(west.p.Dgen_df$LinearFst~west.p.Dgeo_df$km),col="red")

west.p.ibdL <- mantel.randtest(west.p.Dgen,west.p.DgeoL)
west.p.ibdL
par(mfrow=c(1,2),mar=c(5,5,5,1)+.1)
plot(west.p.ibdL)
plot(west.p.DgeoL_df$kmL,west.p.Dgen_df$LinearFst,xlab="ln(km)",ylab="Linearized FST",main="west only")
abline(lm(west.p.Dgen_df$LinearFst~west.p.DgeoL_df$kmL),col="red")

west.p.Dgen_melt <- na.omit(melt(west.p.Dgen_wide[,-1]))
west.p.Dgeo_melt <- na.omit(melt(west.p.Dgeo_wide[,-1]))
west.p.DgeoL_melt <- na.omit(melt(west.p.DgeoL_wide[,-1]))

west.p.lm_ibd <- lm(west.p.Dgen_melt$value~west.p.Dgeo_melt$value)
summary(west.p.lm_ibd)
west.p.lm_ibdL <- lm(west.p.Dgen_melt$value~west.p.DgeoL_melt$value)
summary(west.p.lm_ibdL)


