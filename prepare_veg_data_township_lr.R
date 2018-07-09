#--------------------------------------------------------------------------------------------------------------------------------------------
#In this script vegetation is prepared to estimate the spatial model
#--------------------------------------------------------------------------------------------------------------------------------------------
library(fields)
library(rstan)
library(stepps)

#-------------------------------------------------------------------------------------------------------------------
setwd('~/workflow_stepps_calibration/vegetation/')
help.fun.loc <- 'helper_funs/'
data.loc <- 'data_township/'
plot.loc <- 'plots/'
#-------------------------------------------------------------------------------------------------------------------------------------------
#load data that was used for calibration (is the same vegetation data)
#-------------------------------------------------------------------------------------------------------------------------------------------
veg_coords <- readRDS(paste(data.loc,'vegetation_coords_lr.RDS',sep=''))



#----------------------------------------------------------------------------------------------------------------------
#load township count data
township_data <- readRDS('data_township/tree_counts.RDS')
township_coords <-township_data[,c("Longitude","Latitude")] 
TS_matrix <- readRDS('data_township/township_weights_lr.RDS')
y <- township_data[,-c(1:7)]
N_township <- nrow(y)

#--------------------------------------------------------------------------------------------------------------------------------------------------
# extract coordinates of knots through k-means 
# use k-means because it estimates a center of the coordinates
#-----------------------------------------------------------------------------------------------------------------------------------------------
clust_n <- 120
knot_coords = kmeans(veg_coords, clust_n)$centers 
knot_coords = unname(knot_coords)
N_knots     = dim(knot_coords)[1]
K <- ncol(y)
#d <- rdist(veg_coords,veg_coords)/10^6

#-------------------------------------------------------------------------------------------------------------------
plot(veg_coords,pch = 16,cex = 0.5)
points(knot_coords,col=2,pch =16)
#-------------------------------------------------------------------------------------------------------------------------------
# new d_inter
#-------------------------------------------------------------------------------------------------------------------------------
d_inter <- rdist(veg_coords,knot_coords)/10^6
d_knots <- rdist(knot_coords,knot_coords)/10^6
##############################################################################################################################
N = dim(veg_coords)[1]
x = matrix(1, nrow=N, ncol=1)
N_p = N

temp = qr(x)
Q = qr.Q(temp)
R = qr.R(temp)

P = Q %*% t(Q)
# M = diag(N) - P

if (all(P-P[1,1]<1.0e-12)){
  P = P[1,1]
  N_p = 1
}
##########################################################################################################################
## chunk: save data to file
##########################################################################################################################
taxa <- colnames(y)
veg_coords_agg <- veg_coords

stan_rdump(c('K', 'N', 'N_knots', 
             'y','N_township','TS_matrix', 
             'd_knots', 'd_inter',
             'P', 'N_p'), 
           file=paste(data.loc,'/township_data_', K, '_taxa_', N, '_cells_', N_knots, '_knots.dump',sep=""))

# K = number of species, N = number of vegetation grid points,N_knots = number of knots
# y = vegetation data (rounded to real number), d_knots = distance matrix among knots, d_inter = distance matrix between knots and vegetation grid cells 
# N_p = N, P = N x N matrix all values the same 

save(K, N, N_knots, 
     y, N_township,TS_matrix, 
      d_knots, d_inter, 
     P, N_p,
     veg_coords_agg, knot_coords, taxa,
     file=paste(data.loc,'/township_data_', K, '_taxa_', N, '_cells_', N_knots, '_knots.rdata',sep=""))

