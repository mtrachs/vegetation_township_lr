#--------------------------------------------------------------------------------------------------------------------
#this is actually the code that does something useful
library(rioja)
library(maptools)
library(rgdal)
library(stepps)
dataDir <- '~/github_changed_files/composition/data/'
load(file = file.path(dataDir, paste0('intersection_eastern_twonships_lr.Rda')))


#------------------------------------------------------------------------------------------------------------
#load coordinates of data used 
load('~/workflow_stepps_calibration/calibration/data/elicitation_neus_certainty_median_79_sites_only_abies_new_species.RData')
east <- sort(unique(veg_coords[,'meters.east']))
north <- sort(unique(veg_coords[,'meters.north']))


#now we have to go to low resolution taking into account coordiantes used for lr data
#x.coords
#y.coords

veg_coords_agg.x <- seq(min(east) +8000,max(east),24000)
veg_coords_agg.y <- seq(min(north) +8000,max(north),24000)
veg_coords_agg <- expand.grid(veg_coords_agg.x,veg_coords_agg.y)

dist.matrix <- paldist2(veg_coords_agg,veg_coords,dist.method = 'euclidean')
coord.index.use <- unlist(apply(dist.matrix,2,function(z) which(z==0)))
veg_coords_agg <- veg_coords_agg[coord.index.use,] # is good confirmed by map!!

dist.veg.veg.agg <- paldist2(veg_coords_agg,veg_coords,dist.method = 'euclidean')
min.index <- apply(dist.veg.veg.agg,2,function(y) which.min(y))
r_agg <- stats::aggregate(r,list(cell=min.index),FUN= mean)

#veg_coords
#-------------------------------------------------------------------------------------------------------------
us.shp <- readShapeLines('~/workflow_stepps_calibration/prediction/data/map_data/us_alb.shp',proj4string=CRS('+init=epsg:3175'))


#------------------------------------------------------------------------------------------------------------
#take the initial poly stack data, dimension is x-coordinates, y-coordinates, n_townships 
dim.pstack <- dim(poly.stack)
#transform this data to a dimension of x-coordinates * y-coorindates, n_towhnsips, then we are sure we have the right coordinates
data.ne <- matrix(nrow = dim.pstack[1]*dim.pstack[2],ncol = dim.pstack[3])  

for(i in 1:nTowns) {
  stack <- poly.stack[[i]]@data@values
  data.ne[,i] <- stack/sum(stack)
}

coordinates.tot <- coordinates(poly.stack)

#plot(coordinates.tot)
#points(coordinates.tot[data.ne[,1]>0,],col='orange',pch = 16)


# in fact intersect total coordinates with coordinates available in
rows.use1 <- matrix(ncol=1,nrow=nrow(veg_coords_agg))

for (xxx in 1:nrow(veg_coords_agg)){
  rows.use1[xxx] <- which((coordinates.tot[,1]%in% veg_coords_agg[xxx,1]) & (coordinates.tot[,2]%in%veg_coords_agg[xxx,2]))
}

# dist.coords <- paldist2(veg_coords,coordinates.tot,dist.method='euclidean')
# 
# rows.use <- apply(dist.coords,2,function(x) sum(x==0))
# 
# sum(rows.use>0)#6796 the correct number


coordinates.neus <- coordinates.tot[rows.use1,] #coordinates.tot[rows.use>0,]
data.ne.use <- data.ne[rows.use1,]#data.ne[rows.use>0,]


plot(us.shp,xlim=range(east),ylim=range(north))
for(i in 1:nTowns)
  points(coordinates.neus[data.ne.use[,i]>0,],col='orange',pch = 16,cex = 0.25)

for(i in which(colSums(data.ne.use)<1))
  points(coordinates.neus[data.ne.use[,i]>0,],col='blue',pch = 16,cex = 0.25)

#should not plot those because they are not in the domain
for(i in which(colSums(data.ne.use)==0))
  points(coordinates.neus[data.ne.use[,i]>0,],col=2,pch = 16,cex = 0.25)

#--------------------------------------------------------------------------------------------------------------------
#have to find those that are not in the domain 
#first have to load counts

#load township survey data
#should move this data to other folder, not good to pull dat from different folders!!!
township_data <- read.csv('~/vegetation_data/SetTreeComp_Northeast_Level1_v1.0.csv',head=TRUE)


coordinates.count.data <- township_data[,c("Longitude","Latitude")] 

sputm <- SpatialPoints(coordinates.count.data, proj4string=CRS("+init=epsg:4326"))
spgeo <- spTransform(sputm, CRS("+init=epsg:3175"))
coord_cd_us <- spgeo@coords


dist_km <- matrix(nrow=nrow(coord_cd_us),ncol=1)

for(i in 1:nrow(coord_cd_us)) {
  dist_km1 <- paldist2(matrix(coord_cd_us[i,],ncol=2),matrix(ncol=2,coordinates.neus[data.ne.use[,i]>0,]),dist.method='euclidean')/1000
  dist_km[i] <- min(dist_km1)
}

# use data in domain
tdd <- township_data[colSums(data.ne.use)!=0,]
weight.neus.d <- data.ne.use[,colSums(data.ne.use)!=0] 
weight.neus.d <- apply(weight.neus.d,2,function(x) x/sum(x))
#tcdd <- tdd[,-c(1:7)] # count data

#aggregate counts into correct units
tree_table <- readr::read_csv('~/github_changed_files/composition/data/veg_trans_edited.csv')
#for some reason this code does also reshuffel coordinates
tree_trans <- translate_taxa(tdd, tree_table,id_cols = colnames(tdd)[1:7])
#have to reorder data to match 
tree_trans <- tree_trans[order(tree_trans$ID),]

saveRDS(tree_trans,file='~/workflow_stepps_calibration/vegetation/data_township/tree_counts.RDS')
saveRDS(weight.neus.d,file='~/workflow_stepps_calibration/vegetation/data_township/township_weights_lr.RDS')
saveRDS(r_agg,file='~/workflow_stepps_calibration/vegetation/data_township/vegetation_lr.RDS')
saveRDS(coordinates.neus,file='~/workflow_stepps_calibration/vegetation/data_township/vegetation_coords_lr.RDS')
