#' R function to perform environmental systematic sampling 
#'@description  R function to perform environmental systematic sampling 
#'@usage eSample(envData,nExpect,plot=TRUE,saveShape=TRUE,nf,lowerLim,upperLim)
#'@param envData  Either a SpatialPointsDataFrame, a SpatilPixelsDataFrame, or a SpatialGridDataFrame as defined in package sp or a RasterStack or RasterBrick as defined in package raster. Note that this object can contain a mixture of variables type (quantitative, factor and ordered)
#'@param nExpect Numeric number of desired sampling points. Note that sometime the exact expected number can not be reached because this function tries to find a grid that best matches the expected number
#'@param plot A logic indicating whether or not you would like to have a graphical representation of the environmental systematic sampling (default=TRUE)
#'@param saveShape A logic indcating whether or not you want to save the geographical coordinates of the pixels corresponding to this environmental systematic sampling as a shapefile format (default=TRUE). The generated shapefile will be saved in your working directory
#'@param nf A numeric value indicating the number of ordination axes retained in the analysis, this function manages currently only 2 or 3 dimentions
#'@param lowerLim Numeric value of probability in [0,1] that can be used to produce sample quantiles corresponding to the given lower limit probability. This may be used to reduce the effect of extremes observations
#'@param upperLim Numeric value of probability in [0,1] that can be used to produce sample quantiles corresponding to the given upper limit probability. This may be used to reduce the effect of extremes observations
#'@return Returns a list containing 4 objects: GeoSamples (The geographical coordinates of the pixels corresponding to the environmental systematic sampling ); EnvSamples(The coordinates in the reduced environmental space of pixels corresponding to the environmental systematic sampling ); EnvGrid (The coordinates of the grid in the reduced environmental space)
#'@details The environmental systematic sampling consists in reducing the environmental space into 2 or 3 dimensions using an ordination method (the dudi.mixt method in ade4 package) as a first step. Thereafter convex hull will be created around the reduced environmental space. Then, a multidimensional grid will be created inside the convex hull. The obtained grid represents the perfect configuration required to adequately survey the environmental space. The last step is subsequently to seek the closest pixel to this ideal configuration. This is achieved by searching the nearest neighbour between each grid point and each pixel in the environmental space based on Euclidean distances
#'@references Tarek Hattab, Michael Ewald, Sandra Skowronek, Raf Aerts, Carol X. Garzón-López, Hélène Horen, Boris Brasseur, Emilie Gallet-Moron, Fabien Spicher, Guillaume Decocq, Hannes Feilhauer, Olivier Honnay, Pieter Kempeneers, Sebastian Schmidtlein, Ben Somers, Ruben Van De Kerchove, Duccio Rocchini and Jonathan Lenoir (In prep). An unified framework to model the potential and realized distributions of invasive species within the invaded range

#'@examples
#'library(raster)
#'envData<-getData('worldclim', var='bio', res=10)
#'envData<-crop(envData,extent(-10,45,20,75))
#'plot(envData)
#'par(mfrow=c(1,2))
#'Mysampling1<-eSample(envData,nExpect=50,plot=TRUE,saveShape=TRUE,nf=2,lowerLim=0.00001,upperLim=0.99999)
#'plot(envData[[1]])
#'plot(Mysampling1[[1]],add=TRUE,col=2,pch=19)
#'par(mfrow=c(1,2))
#'Mysampling2<-eSample(envData,nExpect=50,plot=TRUE,saveShape=TRUE,nf=2,lowerLim=0.1,upperLim=0.9)
#'plot(envData[[1]])
#'plot(Mysampling2[[1]],add=TRUE,col=2,pch=19)
#'Mysampling3<-eSample(envData,nExpect=50,plot=TRUE,saveShape=TRUE,nf=3,lowerLim=0.001,upperLim=0.999)
#'plot(envData[[1]])
#'plot(Mysampling3[[1]],add=TRUE,col=2)
#'@export

eSample<-function(envData,nExpect=NULL,plot=TRUE,saveShape=TRUE,nf=3,lowerLim=0.00005,upperLim=0.99995){
  require(raster,quietly = TRUE)
  require(sp,quietly = TRUE)
  require(ade4,quietly = TRUE)
  require(geometry,quietly = TRUE)
  require(pdist,quietly = TRUE)
  require(rgdal,quietly = TRUE)
  require(rgl,quietly = TRUE)
  
  if(is.null(nExpect)){stop("you must provide the nExpect argument")}
  if(!(nf %in% c(2,3))) {stop("nf must be equal to either 2 or 3") }
  if (lowerLim >1 | lowerLim < 0) {stop("lowerLim must be between 0 and 1") }
  if (upperLim >1 | upperLim < 0) {stop("upperLim must be between 0 and 1") }
  if (upperLim < lowerLim) {stop("lowerLim must be smaller than upperLim") }
  if(class(envData) %in% c("SpatialGridDataFrame")){x<-envData@data
                                                    rownames(x)<-seq(1:dim(x)[1])
                                                    XY<-coordinates(envData)
                                                    rownames(XY)<-rownames(x)
                                                    x<-na.omit(x)
                                                    XY<-XY[rownames(x),]}
  if(class(envData) %in% c("SpatialPointsDataFrame","SpatilPixelsDataFrame")){x<-envData@data
                                                                              rownames(x)<-seq(1:dim(x)[1])
                                                                              XY<-coordinates(envData)}
  if(class(envData) %in% c("RasterStack","RasterBrick")){x<-values(envData)
                                                         XY<-coordinates(envData)
                                                         rownames(x)<-seq(1:dim(x)[1])
                                                         rownames(XY)<-rownames(x)
                                                         x<-na.omit(x)
                                                         XY<-XY[rownames(x),]}
  if(!(class(envData) %in% c("SpatialPointsDataFrame","SpatialGridDataFrame","RasterStack","RasterBrick"))){stop("The envData object must be either a SpatialPointsDataFrame or a SpatialGridDataFrame or a RasterStack or a RasterBrick" )}
  
  if (nExpect > dim(x)[1]) {stop("The number of sampling points required is greater than the number of pixels in the environmental data")}
  
  Ordi<-dudi.mix(x,nf=nf,scannf=F)
  id<-rownames(x)
  x<-Ordi$li
  rownames(x)<-id
  lim.lower<-apply(x,2, function(x) quantile(x,probs=lowerLim))
  lim.upper<-apply(x,2, function(x) quantile(x,probs=upperLim))
  for(i in 1:dim(x)[2]) {x<-x[which(x[,i]>lim.lower[i] & x[,i]<lim.upper[i]),]}
  XY<-XY[rownames(x),]
  conv<-t(convhulln(x))
  hull<-cbind(x=x[conv,1],y=x[conv,2],z=x[conv,3])
  bb<-cbind(min=apply(hull,2,min),max=apply(hull,2,max))
  
  inhull <- function(testpts, calpts, hull=convhulln(calpts), tol=mean(mean(abs(calpts)))*sqrt(.Machine$double.eps)) {
    require(geometry,quietly = TRUE)
    require(MASS,quietly = TRUE)
    calpts <- as.matrix(calpts)
    testpts <- as.matrix(testpts)
    p <- dim(calpts)[2]
    cx <- dim(testpts)[1]
    nt <- dim(hull)[1]
    nrmls <- matrix(NA, nt, p)
    degenflag <- matrix(TRUE, nt, 1)
    for (i in 1:nt) {nullsp <- t(Null(t(calpts[hull[i,-1],] - matrix(calpts[hull[i,1],],p-1,p, byrow=TRUE))))
                     if (dim(nullsp)[1] == 1) { nrmls[i,] <- nullsp
                                                degenflag[i] <- FALSE}}
    if(length(degenflag[degenflag]) > 0) warning(length(degenflag[degenflag])," degenerate faces in convex hull")
    nrmls <- nrmls[!degenflag,]
    nt <- dim(nrmls)[1]
    center = apply(calpts, 2, mean)
    a <- calpts[hull[!degenflag,1],]
    nrmls <- nrmls/matrix(apply(nrmls, 1, function(x) sqrt(sum(x^2))), nt, p)
    dp <- sign(apply((matrix(center, nt, p, byrow=TRUE)-a) * nrmls, 1, sum))
    nrmls <- nrmls*matrix(dp, nt, p)
    aN <- diag(a %*% t(nrmls))
    val <- apply(testpts %*% t(nrmls) - matrix(aN, cx, nt, byrow=TRUE), 1,min)
    val[abs(val) < tol] <- 0
    as.integer(sign(val))}
  
  createGrid<-function(nExpect=nExpect,nsig=2,bb,x=x){
    n = nExpect; offset = rep(0.1,  nrow(bb)); nProd<-0
    pb <- txtProgressBar(min=0,max=nExpect,style = 3)
    while(nProd<nExpect){pw = 1/nrow(bb)
                         cellsize = signif((prod(apply(bb, 1, diff))/n)^pw, nsig)
                         cellsize = min(cellsize, min(apply(bb, 1, diff)))
                         cellsize<-rep(cellsize,nf)
                         min.coords = pmax(bb[, 1], signif(bb[, 1] + offset * cellsize, nsig))
                         expand.grid.arglist = list()
                         for (i in 1:nrow(bb)) {name = paste("x", i, sep = "")
                                                sign = ifelse(min.coords[i] < bb[i, 2], 1, -1)
                                                expand.grid.arglist[[name]] = seq(min.coords[i], bb[i,2], sign * cellsize[i])}
                         xyz = do.call(expand.grid, expand.grid.arglist)
                         selection<-inhull(testpts=xyz, calpts=x, hull=convhulln(x))
                         xyz<-xyz[which(selection==1),]
                         distances<-pdist( xyz,x)
                         distances<-as.matrix(distances)
                         XY.samples<-NULL
                         for (i in 1:dim(xyz)[1]) {XY.samples<-rbind (XY.samples,XY[rownames(x)[which.min(distances[i,1:dim(x)[1]])],])}
                         EnvSamples<-NULL
                         for (i in 1:dim(xyz)[1]) {EnvSamples<-rbind (EnvSamples,x[rownames(x)[which.min(distances[i,1:dim(x)[1]])],])}
                         XY.samples<-as.data.frame(XY.samples)
                         XY.samples$id<-seq(1:dim(XY.samples)[1])
                         coordinates(XY.samples)<-~x+y
                         XY.samples<-remove.duplicates(XY.samples)
                         nProd<-dim(XY.samples)[1]
                         n<-n+1
                         Sys.sleep(0.1)
                         setTxtProgressBar(pb,nProd )}
    close(pb)
    output<-list(GeoSamples=XY.samples,EnvSamples=EnvSamples,EnvGrid=xyz)
    return(output)
  }
  print(paste("-----------------Optimization of the",nf, "D grid size---------------------"))
  output<-createGrid(nExpect=nExpect,nsig=2,bb=bb,x=x)
  print(paste("-----------------",dim(output[[1]])[1]," points found---------------------",sep=""))
  if(saveShape==T){writeOGR(output[[1]], getwd(), "EnvSysSample", driver="ESRI Shapefile",overwrite_layer = T)}
  if(plot==T & nf==3){rgl::open3d()
                      rgl::plot3d(x[,1], x[,2], x[,3], r = 0.2,xlab="Axis 1",ylab="Axis 2",zlab="Axis 3")
                      rgl::plot3d(output[[2]][,1],output[[2]][,2],output[[2]][,3], col=2, type ="s",r=.1,add=T)
                      rgl::rgl.triangles(x[conv,1],x[conv,2],x[conv,3],col=2,alpha=0.2)
                      moncube<-rgl::cube3d(col = "red")
                      moncube$vb[1,] <- moncube$vb[1,]*output[[4]][1]
                      moncube$vb[2,] <- moncube$vb[2,]*output[[4]][2]
                      moncube$vb[3,] <- moncube$vb[3,]*output[[4]][3]
                      for(i in 1:dim(output[[3]])[1]){rgl::wire3d( rgl::translate3d( moncube, output[[3]][i,1], output[[3]][i,2],output[[3]][i,3]))}}
 
  if(plot==T & nf==2){  plot(x[,1],x[,2],cex=0.3,xlab="Axis 1",ylab="Axis 2")
                        points(output[[2]][,1],output[[2]][,2],col=2,cex=0.7,pch=16)
                        x2<-x[,c(1,2)]
                        coordinates(x2)<-~Axis1+Axis2
                        ch <- chull(coordinates(x2))
                        coords <- coordinates(x2)[c(ch, ch[1]), ]
                        sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
                        plot(sp_poly,add=T,bg=2,border=2)   }
  return(output)
} # end of function eSample


