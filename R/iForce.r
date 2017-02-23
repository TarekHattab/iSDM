#' R function to compute a negative exponential dispersal kernel
#'@description  R function to compute a negative exponential dispersal kernel 
#'@usage  iForce(occData,coords,a,envData,binary=TRUE,longlat=NULL)
#'@param occData  Either a SpatialPointsDataFrame as defined in package sp, a data.frame or a matrix object containing species data 
#'@param coords Optional 2 columns matrix containing the X and Y coordinates of occData (only consider if occData is a a data.frame or a matrix object)
#'@param a An integer between 0 and 1 that controls the form of the dispersal kernel 
#'@param envData Either a SpatialPixelsDataFrame or a SpatialGridDataFrame as defined in package sp or a Raster or RasterLayer as defined in package raster. This object will be used to determine the area for which the negative exponential dispersal kernel will be calculated
#'@param binary A logic indicating whether the occData object contains presences/absences data or only the geographical coordinates of presences data (default=TRUE)
#'@param longlat A logic indicating whether point coordinates are projected (longitude-latitude decimal degrees) or not (for plane) (default=TRUE)
#'@return Returns a RasterLayer object of the negative exponential dispersal kernel 
#'@references Tarek Hattab, Carol Ximena Garzon Lopez, Michael Ewald, Sandra Skowronek, Raf Aerts, Helene Horen, Boris Brasseur, Emilie Gallet-Moron, Fabien Spicher, Guillaume Decocq, Hannes Feilhauer, Olivier Honnay, Pieter Kempeneers, Sebastian Schmidtlein, Ben Somers, Ruben Van De Kerchove, Duccio Rocchini and Jonathan Lenoir (Accpeted). A unified framework to model the potential and realized distributions of invasive species within the invaded range. Diversity and Distributions.
#'@export
#'@examples
#'library(raster)
#'library(sp)
#'envData<-getData('worldclim', var='bio', res=10)
#'envData<-crop(envData,extent(-8,15,38,55))
#'envData<-envData[[1]]
#'
#' #Generate randomly a SpatilaPointsDataFrame containing occurrences
#'xy<-coordinates(envData)[sample(which(is.na(values(envData))==FALSE),100),]
#'occ<-ifelse(xy[,2]>50,0,1)
#'occData<-SpatialPointsDataFrame(coords=xy,data=as.data.frame(occ), 
#'proj4string = CRS(proj4string(envData)))
#'
#'par(mfrow=c(3,3),mar=c(1.5,1.5,1.5,1.5))
#'for (a in c(0.01,0.03,0.05,0.06,0.09,0.1,0.3,0.6,0.9)){
#'propagule<-iForce(occData,envData=envData,a=a,binary=TRUE,longlat=TRUE)
#'plot(propagule,main=paste("a = ",a))
#'plot(occData,col=ifelse(occData@data[,1]==1,1,0),add=TRUE,cex=0.3)}

iForce<-function(occData,coords=NULL,a=NULL,envData,binary=TRUE,longlat=NULL){

  if(is.null(longlat)){stop("you must provide a longlat argument (TRUE/FALSE)")}

  if(is.null(a)){stop("you must provide the a argument")}

  if(!(class(occData) %in% c("data.frame","matrix","SpatialPointsDataFrame"))) {stop("The occData object must be either a SpatialPointsDataFrame , a data.frame or a matrix object")}

  if(dim(occData)[2]>1){stop("The number of columns in occData is greater than 1")}

  if(class(occData)!="SpatialPointsDataFrame" & is.null(coords)==T) {stop("the coords object containing the X and Y coordinates of occData must   be provided")}

  if( is.null(coords)==F & class(occData) %in% c("data.frame","matrix") ){ if (dim(coords)[1]!= dim(occData)[1]) {stop("coords and occData have a  different size ")}}

  if( is.null(coords)==F & class(occData) %in% c("data.frame","matrix")){occData<-as.data.frame(cbind(coords,occData))
                                                                          names(occData)<-c("x","y","SP")
                                                                          sp::coordinates(occData)<-~x+y}

  if(!all(is.element(unique(occData@data[,1]), c(0,1)))& binary==T) {stop("occData contains values different from 0 and 1")}

  if(class(occData)=="SpatialPointsDataFrame" & binary==T) {presence<-sp::coordinates(occData[which(occData@data[,1]==1),])}

  if(class(occData)=="SpatialPointsDataFrame" & binary==F) {presence<-sp::coordinates(occData)}

  if(!(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame","Raster","RasterLayer"))){stop("The envData  object must be either a SpatialPixelsDataFrame a SpatialGridDataFrame a RasterStack or a RasterBrick" )}

  if(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame")) {envData<-raster::raster(envData)}
   
  XY<-sp::coordinates(envData)
  distance<-raster::pointDistance(presence,XY,lonlat=longlat)
  a<-max(distance)*a
  dist.sq<-exp(-distance/a)
  propagule<-colSums(dist.sq)
  propagule<-raster::rasterFromXYZ(cbind(XY,propagule),res=raster::res(envData),digits=1e-20)
  propagule<-raster::mask(propagule,envData)

  return(propagule)
}# end of function iForce
