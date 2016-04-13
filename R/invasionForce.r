#' R function to compute a negative exponential dispersal kernel
#'@description  R function to compute a negative exponential dispersal kernel 
#'@usage  invasionForce(occData,coords,a,envData,binary=TRUE,longlat=NULL)
#'@param occData May be a SpatialPointsDataFrame as defined in package sp, a data.frame or  a matrix object containing species data 
#'@param coords Optional 2 columns matrix containing the X and Y coordinates of occData(only consider if occData is a a data.frame or  a matrix object)
#'@param a An integer between 0 and 1 that controls the form of the dispersal kernel 
#'@param envData May be a SpatialPointsDataFrame or a SpatialGridDataFrame as defined in package sp or a RasterStack or RasterBrick as defined in package raster. This object will be used to determine the area for which the negative exponential dispersal kernel will be calculated
#'@param binary logical TRUE if the occData object contains presences/absences data and FALSE if occData contains only the geographical coordinates of presences data
#'@param longlat Logical If TRUE, coordinates should be in degrees; else they should represent planar (Euclidean) space (e.g. units of meters)
#'@return Returns a RasterLayer object of the negative exponential dispersal kernel 
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
#'occData<-SpatialPointsDataFrame(coords=xy,data=as.data.frame(occ), proj4string = CRS(proj4string(envData)))
#'
#'par(mfrow=c(3,3),mar=c(1.5,1.5,1.5,1.5))
#'for (a in c(0.01,0.03,0.05,0.06,0.09,0.1,0.3,0.6,0.9)){
#'propagule<-invasionForce(occData,envData=envData,a=a,binary=TRUE,longlat=TRUE)
#'plot(propagule,main=paste("a = ",a))
#'plot(occData,col=ifelse(occData@data[,1]==1,1,0),add=TRUE,cex=0.3)}

invasionForce<-function(occData,coords=NULL,a=NULL,envData,binary=TRUE,longlat=NULL){

  require(raster,quietly = TRUE)
  require(sp,quietly = TRUE)

  if(is.null(longlat)){stop("you must provide a longlat argument (TRUE/FALSE)")}

  if(is.null(a)){stop("you must provide the a argument")}

  if(!(class(occData) %in% c("data.frame","matrix","SpatialPointsDataFrame"))) {stop("The occData object must be either a SpatialPointsDataFrame , a data.frame or a matrix object")}

  if(dim(occData)[2]>1){stop("The number of columns in occData is greater than 1")}

  if(class(occData)!="SpatialPointsDataFrame" & is.null(coords)==T) {stop("the coords object containing the X and Y coordinates of occData must   be provided")}

  if( is.null(coords)==F & class(occData) %in% c("data.frame","matrix") ){ if (dim(coords)[1]!= dim(occData)[1]) {stop("coords and occData have a  different size ")}}

  if( is.null(coords)==F & class(occData) %in% c("data.frame","matrix")){occData<-as.data.frame(cbind(coords,occData))
                                                                          names(occData)<-c("x","y","SP")
                                                                          coordinates(occData)<-~x+y}

  if(!all(is.element(unique(occData@data[,1]), c(0,1)))& binary==T) {stop("occData contains values different from 0 and 1")}

  if(class(occData)=="SpatialPointsDataFrame" & binary==T) {presence<-coordinates(occData[which(occData@data[,1]==1),])}

  if(class(occData)=="SpatialPointsDataFrame" & binary==F) {presence<-coordinates(occData)}

  if(!(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame","Raster","RasterLayer"))){stop("The envData  object must be either a SpatialPointsDataFrame a SpatialGridDataFrame a RasterStack or a RasterBrick" )}

  if(class(envData) %in% c("SpatialGridDataFrame","SpatialPixelsDataFrame")) {envData<-raster(envData)}
   
  XY<-coordinates(envData)
  distance<-pointDistance(presence,XY,lonlat=longlat)
  a<-max(distance)*a
  dist.sq<-exp(-distance/a)
  propagule<-colSums(dist.sq)
  propagule<-rasterFromXYZ(cbind(XY,propagule),res=res(envData),digits=1e-20)
  propagule<-mask(propagule,envData)

  return(propagule)
}# end of function invasionForce
