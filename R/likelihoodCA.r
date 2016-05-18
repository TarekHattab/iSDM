#' R function to compute the likelihood of detecting contingent absences within a network of presence-absence data
#'@description R function to compute the likelihood of detecting contingent absences
#'@usage likelihoodCA(occData,coords=NULL,envData, longlat=TRUE,nf=5)
#'@param occData May be a SpatialPointsDataFrame as defined in package sp, a data.frame or  a matrix object containing species data (a single species) in binary format (ones for presences, zeros for absences) 
#'@param coords Optional 2 columns matrix containing the X and Y coordinates of occData (only consider if occData is a a data.frame or  a matrix object) 
#'@param envData An abject containing your explanatory variables. May be a SpatialPointsDataFrame or SpatialGridDataFrame as defined in package sp or a RasterStack or RasterBrick as defined in package raster. Note that this object can contain a mixed type variables (quantitative, factor and ordered).
#'@param longlat Logical TRUE if point coordinates are not projected(longitude-latitude decimal degrees)(+proj=longlat) and FALSE for plane (lonlat=FALSE)
#'@param nf If envData contains any factor and ordered variables, the Hill and Smith method will be used to perform an ordination of explanatory variables. nf corresponds in this case to the number of kept axes for the caluluation of mahalanobis distances. 
#'@details This function campute the likelihood of detecting contingent absences within a network of presence-absence data. It's based on the idea that absences data that are environmentally close but geographically distant to presences data are more likely to be contingent absences. 
#'@return Returns a SpatialPointsDataFrame containing the likelihood values  
#'@export
#'@references Tarek Hattab, Ruben Van De Kerchove, Ben Somers, Boris Brasseur, Carol Ximena Garz?n L?pez, Duccio Rocchini, Emilie Gallet-Moron, Fabien Spicher, Guillaume Decocq, Hannes Feilhauer, Hélène Horen, Jens Warrie, Michael Ewald, Olivier Honnay, Pieter Kempeneers, Raf Aerts, Sandra Skowronek, Sebastian Schmidtlein and Jonathan Lenoir (In prep).An unified framework to model the potential and realized distributions of invasive species within the invaded range
#'@examples
#'library(raster)
#'library(sp)
#'library(virtualspecies)
#'library(colorRamps)
#'
#'envData<-getData('worldclim', var='bio', res=10)
#'envData<-crop(envData,extent(-8,15,38,55))
#'
#' #Generate virtual species responses with formatfunctions 
#'my.parameters <- formatFunctions(bio1 = c(fun = "dnorm", mean = 140, sd = 40), bio5 = c(fun = "dnorm",
#'mean = 230, sd = 70),bio6 = c(fun = "dnorm",mean = 10, sd = 40))
#'
#' #Generate a virtual species distributions with responses to environmental variables
#'potential.dist <- generateSpFromFun(envData[[c(1,5,6)]], my.parameters)
#'
#' #Limit a virtual species distribution to a defined area. 
#' #It will thus generate a species which is not at the equilibrium with its environment
#'realized.dist<-limitDistribution(x=potential.dist$suitab.raster, area=extent(-8,15,38,48))
#'
#' #Generate a random presence absence dataset from the realized distribution 
#'# using a probability threshold of 0.5
#'Presence<-coordinates(realized.dist$occupied.area
#')[sample(which(values(realized.dist$occupied.area)>0.5),300),]
#'Absence<-coordinates(realized.dist$occupied.area
#')[sample(which(values(realized.dist$occupied.area)<0.5),300),]
#'occData<-as.data.frame(rbind(cbind(Presence,SP=rep(1,300)),cbind(Absence,SP=rep(0,300))))
#'coordinates(occData)<-~x+y
#'proj4string(occData)<-proj4string(envData)
#'
#' # Calculat  the likelihood of detecting contingent absences
#'likelihood<-likelihoodCA(occData=occData,envData=envData[[c(1,5,6)]],longlat=TRUE)
#'
#' # Display the results
#'par(mfrow=c(1,2),mar=c(2,2.5,2,2.5))
#'plot(realized.dist$occupied.area,main="Realized distribution")
#'plot(occData,col=ifelse(occData$SP==1,2,1),add=TRUE,pch=19,cex=0.8)
#'plot(potential.dist$suitab.raster,main="Potential distribution")
#'
#'scatterCol<-function(x){
#'x<-(x-min(x))/(max(x)-min(x))
#'colorFunction <- colorRamp(matlab.like(100))
#'zMatrix <- colorFunction(x)
#'zColors <- rgb(zMatrix[,1], zMatrix[,2], zMatrix[,3], maxColorValue=255)
#'return(zColors)}
#'
#'points(likelihood,pch=21, col=1,bg=scatterCol(likelihood@data[,"LCA"]),cex=1)


likelihoodCA<-function(occData,coords=NULL, envData, longlat=TRUE,nf=5){
  
  require(raster,quietly = TRUE)
  require(sp,quietly = TRUE)
  require(ade4,quietly = TRUE)
  require(geosphere,quietly = TRUE)
  require(maptools,quietly = TRUE)
  require(pdist,quietly=TRUE)
  require(ade4,quietly=TRUE)
  
  if(is.null(longlat)){stop("you must provide a longlat argument (TRUE/FALSE)")}
  
  if (!(class(occData) %in% c("data.frame","matrix","SpatialPointsDataFrame"))) {stop("The occData object must be either a SpatialPointsDataFrame , a data.frame or a matrix")}
  if (dim(occData)[2]>1){stop("The number of columns in occData is greater than 1")}
  
  if (class(occData)!="SpatialPointsDataFrame" & is.null(coords)==T) {stop("the coords object containing the X and Y coordinates of occData must be provided")}
  if ( is.null(coords)==F & class(occData) %in% c("data.frame","matrix") ){ if (dim(coords)[1]!= dim(occData)[1]) {stop("coords and occData have a different size ")}}
  if ( is.null(coords)==F & class(occData) %in% c("data.frame","matrix")){occData<-as.data.frame(cbind(coords,occData))  
                                                                          names(occData)<-c("x","y","SP")
                                                                          coordinates(occData)<-~x+y} 
  
  if (!all(is.element(unique(occData@data[,1]), c(0,1)))) {stop("occData contains non-binary values")}
  if(class(envData) %in% c("SpatialGridDataFrame,SpatialPointsDataFrame","SpatialPixelsDataFrame")){occEnv<-over(occData,envData)
                                                                                                    occEnv<-cbind(SP=occData@data[,1],occEnv)}
  if(class(envData) %in% c("RasterStack","RasterBrick")){occEnv<-extract(y=occData,x=envData)
                                                         occEnv<-cbind(SP=occData@data[,1],occEnv)}
  if(!(class(envData) %in% c("SpatialPointsDataFrame","SpatialGridDataFrame","SpatialPixelsDataFrame","RasterStack","RasterBrick"))){stop("The envData object must be either a SpatialPointsDataFrame or a SpatialGridDataFrame or a RasterStack or a RasterBrick" )}
  
  presence<-occData[which(occData@data[,1]==1),]
  absence<-occData[which(occData@data[,1]==0),]
  
  if (longlat==TRUE){distGeo<-NULL
                     for (i in 1:length(absence)) {
                       distGeo<- rbind(distGeo,min(distVincentyEllipsoid(presence,absence[i,])))}
                     geoWeight<-(distGeo-min(distGeo))/(max(distGeo)-min(distGeo))}
  
  if(longlat==FALSE){distGeo<-pdist(coordinates(presence),coordinates(absence))
                     distGeo<-as.matrix(distGeo)
                     distGeo<-apply(distGeo,2,FUN=min)
                     geoWeight<-(distGeo-min(distGeo))/(max(distGeo)-min(distGeo))}
  
  absence@data$geoWeight<-geoWeight
  presence@data$geoWeight<-rep(0,dim(presence)[1])
  
  
  if(any(sapply(occEnv,is.factor)[-1])|any(sapply(occEnv,is.ordered)[-1])){ordi<-dudi.mix(occEnv[,-1],scannf=F,nf=nf)
                                                                           occEnv<-cbind(occEnv[,1],ordi$l1)
                                                                           print(paste("The Hill and Smith method was used to perform an ordination of explanatory variables ", nf, " axes were kept for the caluluation of mahalanobis distances",sep=""))}                                                                                       
  mahalanobis.dist<-function (data.x, data.y = NULL, vc = NULL,tol=1e-20) 
  {
    xx <- as.matrix(data.x)
    if (is.null(data.y)) 
      yy <- as.matrix(data.x)
    else yy <- as.matrix(data.y)
    if (is.null(vc)) {
      if (is.null(data.y)) 
        vc <- var(xx)
      else vc <- var(rbind(xx, yy))
    }
    ny <- nrow(yy)
    md <- matrix(0, nrow(xx), ny)
    for (i in 1:ny) {
      md[, i] <- mahalanobis(xx, yy[i, ], cov = vc,tol=tol)
    }
    if (is.null(data.y)) 
      dimnames(md) <- list(rownames(data.x), rownames(data.x))
    else dimnames(md) <- list(rownames(data.x), rownames(data.y))
    sqrt(md)
  }
  
  presenceEnv<-occEnv[which(occEnv[,1]==1),-1]
  absenceEnv<-occEnv[which(occEnv[,1]==0),-1]
  distEnv<-mahalanobis.dist(data.x=presenceEnv,data.y=absenceEnv,tol=1e-20)
  distEnv<-as.matrix(distEnv)
  envWeight<-apply(distEnv,2,FUN=min)
  envWeight<-(1-(envWeight-min(envWeight))/(max(envWeight)-min(envWeight))) 
  
  absence@data$envWeight<-envWeight
  presence@data$envWeight<-rep(0,dim(presence)[1])
  occRes<-spRbind(presence,absence)
  occRes$LCA<-(occRes$envWeight+occRes$geoWeight)/2
  occRes<-occRes[,"LCA"]
  return(occRes)
} # end of function likelihoodCA
