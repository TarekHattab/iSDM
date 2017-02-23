#' R function to compute the probability of detecting dispersal-limited absences within a network of presence-absence data
#'@description R function to compute the probability of detecting dispersal-limited absences
#'@usage pDLA(occData,envData, longlat=TRUE,nf=5,occNative=NULL,envNative=NULL)
#'@param occData A SpatialPointsDataFrame as defined in package sp containing species data (a single species) in binary format (ones for presences, zeros for absences) 
#'@param envData An abject containing your explanatory variables. Either a SpatialPixelsDataFrame or SpatialGridDataFrame as defined in package sp or a RasterStack or RasterBrick as defined in package raster. Note that this object can contain a mixture of variables type (quantitative, factor and ordered).
#'@param longlat A logic indicating whether point coordinates are projected (longitude-latitude decimal degrees) or not (for plane) (default=TRUE)
#'@param nf If envData contains any factor and ordered variables, the Hill and Smith method will be used to perform an ordination of explanatory variables. nf corresponds in this case to the number of kept axes for the caluluation of mahalanobis distances. 
#'@param occNative (optional) a SpatialPoint object containing the occurrences from the native range
#'@param  envNative (optional) an object of the same class as envData containing the explanatory variables from the native range 
#'@details This function campute the probability of detecting dispersal-limited absences within a network of presence-absence data. It's based on the idea that absences data that are environmentally close but geographically distant to presences data are more likely to be dispersal-limited absences. This function allows combining presence data from both the native and the invaded range, note that in this case the set of presence from the native range will only be added to the set of presences from the invaded range when calculating distances between absences and presences in the environnental space as it does not make sens to add presence data from the native range when calculating these distances within the geographical space.
#'@return Returns a SpatialPointsDataFrame containing the probability values  
#'@export
#'@references Tarek Hattab, Carol Ximena Garzon Lopez, Michael Ewald, Sandra Skowronek, Raf Aerts, Helene Horen, Boris Brasseur, Emilie Gallet-Moron, Fabien Spicher, Guillaume Decocq, Hannes Feilhauer, Olivier Honnay, Pieter Kempeneers, Sebastian Schmidtlein, Ben Somers, Ruben Van De Kerchove, Duccio Rocchini and Jonathan Lenoir (Accpeted). A unified framework to model the potential and realized distributions of invasive species within the invaded range. Diversity and Distributions.
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
#'my.parameters <- formatFunctions(bio1 = c(fun = "dnorm", mean = 140, sd = 40),
#' bio5 = c(fun = "dnorm",mean = 230, sd = 70),
#' bio6 = c(fun = "dnorm",mean = 10, sd = 40))
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
#' # Calculat  the probability of detecting dispersal-limited absences
#'probability<-pDLA(occData=occData,envData=envData[[c(1,5,6)]],longlat=TRUE)
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
#'points(probability,pch=21, col=1,bg=scatterCol(probability@data[,"PDLA"]),cex=1)
#'
#'# Example based on occurrences from both the invaded and the native range 
#'envData<-getData('worldclim', var='bio', res=10)
#'envNative<-crop(envData,extent(-90,-70,20,40))
#'envData<-crop(envData,extent(-8,15,38,55))

#'native.dist <- generateSpFromFun(envNative[[c(1,5,6)]], my.parameters)
#'occNative<-as.data.frame(coordinates(native.dist$suitab.raster)
#'[sample(which(values(native.dist$suitab.raster)>0.5),100),])
#'coordinates(occNative)<-~x+y
#'proj4string(occNative)<-proj4string(envData)
#'plot(native.dist,main=" Native range distribution")
#'plot(occNative,add=TRUE,pch=19,cex=0.8)

#'probability<-pDLA(occData=occData,envData=envData[[c(1,5,6)]],longlat=TRUE,
#'occNative=occNative,envNative=envNative[[c(1,5,6)]])

pDLA<-function(occData, envData, longlat=TRUE,nf=5,occNative=NULL,envNative=NULL){
  
  if(is.null(longlat)){stop("you must provide a longlat argument (TRUE/FALSE)")}
  
  if (!(class(occData) %in% c("SpatialPointsDataFrame"))) {stop("The occData object must be a SpatialPointsDataFrame")}
  
  if (dim(occData)[2]>1){stop("The number of columns in occData is greater than 1")}
  
  if (!all(is.element(unique(occData@data[,1]), c(0,1)))) {stop("occData contains non-binary values")}
  
  if (!all(c(is.null(occNative),is.null(envNative)))) {if (class(occNative) != "SpatialPoints") {stop("The occNative object must be a SpatialPoints")}} 
  
  if(!(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame","RasterStack","RasterBrick"))){stop("The envData object must be either a SpatialPixelsDataFrame or a SpatialGridDataFrame or a RasterStack or a RasterBrick" )}
  
  if (!all(c(is.null(occNative),is.null(envNative)))) {if (!(class(envNative)%in% c("SpatialGridDataFrame","SpatialPixelsDataFrame","RasterStack","RasterBrick"))) {stop("The occNative object must be either a SpatialPixelsDataFrame or a SpatialGridDataFrame or a RasterStack or a RasterBrick")}} 
  
  if (!all(c(is.null(occNative),is.null(envNative)))) { if (class(envNative)!=class(envData)) {stop("envNative and envData must belong to the same object class")}}
  
  if (!all(c(is.null(occNative),is.null(envNative)))) {if (any(names(envNative)!=names(envData))) {stop("The variable names in envNative are different from those in envData")}}
  
  if(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame")){occEnv<-sp::over(occData,envData)
                                                                             occEnv<-cbind(SP=occData@data[,1],occEnv)}
  if(class(envData) %in% c("RasterStack","RasterBrick")){occEnv<-raster::extract(y=occData,x=envData)
                                                         occEnv<-cbind(SP=occData@data[,1],occEnv)}
  
  if (!all(c(is.null(occNative),is.null(envNative)))) {if(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame"))
  {occEnv1<-sp::over(occData,envData)
   occEnv1<-cbind(SP=occData@data[,1],occEnv1)
   occEnv2<-sp::over(occNative,envNative)
   occEnv2<-cbind(SP=rep(1,nrow(occEnv2)),occEnv2)
   occEnv<-rbind(occEnv1,occEnv2)}}
  
  if (!all(c(is.null(occNative),is.null(envNative)))) {if(class(envData) %in% c("RasterStack","RasterBrick"))
  {occEnv1<-raster::extract(y=occData,x=envData)
   occEnv1<-cbind(SP=occData@data[,1],occEnv1)
   occEnv2<-raster::extract(y=occNative,x=envNative)
   occEnv2<-cbind(SP=rep(1,nrow(occEnv2)),occEnv2)
   occEnv<-rbind(occEnv1,occEnv2)}}
  
  names.occ<-rownames(occData@data)
  presence<-occData[which(occData@data[,1]==1),]
  absence<-occData[which(occData@data[,1]==0),]
  
  if (longlat==TRUE){distGeo<-NULL
                     for (i in 1:length(absence)) {
                       distGeo<- rbind(distGeo,min(geosphere::distVincentyEllipsoid(presence,absence[i,])))}
                     geoWeight<-(distGeo-min(distGeo))/(max(distGeo)-min(distGeo))}
  
  if(longlat==FALSE){distGeo<-pdist::pdist(sp::coordinates(presence),sp::coordinates(absence))
                     distGeo<-as.matrix(distGeo)
                     distGeo<-apply(distGeo,2,FUN=min)
                     geoWeight<-(distGeo-min(distGeo))/(max(distGeo)-min(distGeo))}
  
  absence@data$geoWeight<-geoWeight
  presence@data$geoWeight<-rep(0,dim(presence)[1])
  if(any(sapply(occEnv,is.factor)[-1])|any(sapply(occEnv,is.ordered)[-1])){
    ordi<-ade4::dudi.mix(occEnv[,-1],scannf=F,nf=nf)
    occEnv<-cbind(occEnv[,1],ordi$l1)
    print(paste("The Hill and Smith method was used to perform an ordination of explanatory variables ",  nf, " axes were kept for the caluluation of mahalanobis distances",sep=""))}    
  
  presenceEnv<-occEnv[which(occEnv[,1]==1),-1]
  absenceEnv<-occEnv[which(occEnv[,1]==0),-1]
  centroide<-colMeans(presenceEnv)
  vc<-stats::var(rbind(absenceEnv,centroide))
  envWeight<-stats::mahalanobis(absenceEnv, centroide, cov = vc,tol=1e-20)
  envWeight<-(1-(envWeight-min(envWeight))/(max(envWeight)-min(envWeight))) 
  absence@data$envWeight<-envWeight
  presence@data$envWeight<-rep(0,dim(presence)[1])
  occRes<-maptools::spRbind(presence,absence)
  occRes$PDLA<-(occRes$envWeight+occRes$geoWeight)/2
  occRes<-occRes[,"PDLA"]
  occRes<-occRes[match(names.occ,rownames(occRes@data)),]
  return(occRes)
  
} # end of function pDLA
