noneucdists = function(masknu,userfn,covariate=NULL) {
  if(!is.element("noneuc",names(covariates(masknu)))) stop("Need a mask covariate called `noneuc'.")
  if(!is.character(covariate)) stop("Argument `covariate' must be a character.")
  if(!is.null(covariate) & !is.element(covariate,names(covariates(masknu)))) 
    stop(paste("No covariate called",covariate,"in the mask"))
  f = match.fun(userfn)
  uNU = range(covariates(masknu)$noneuc)
  loNUind = which(covariates(masknu)$noneuc==uNU[1])[1]
  hiNUind = which(covariates(masknu)$noneuc==uNU[2])[1]
  lonu = covariates(masknu)$noneuc[loNUind]
  hinu = covariates(masknu)$noneuc[hiNUind]
  if(!is.null(covariate)) {
    covcol = which(names(covariates(masknu))==covariate)
    locov = covariates(masknu)[,covcol][loNUind]
    hicov = covariates(masknu)[,covcol][hiNUind]
  }
  # make mask with two points, distance 2 apart
  dets = trap.builder(frame=data.frame(x=c(-2,2),y=c(0,0)),method="all")
  msk = make.mask(dets,buffer=2, spacing=1)
  # calculate dists for lo-lo, lo-hi, hi-lo, hi-hi
  covariates(msk) = data.frame(noneuc=rep(c(lonu,hinu,hinu,lonu,lonu,rep(0,dim(msk)[1]-5))))
  alldists = f(as.matrix(msk),as.matrix(msk),msk)
  dists = c(alldists[1,2],alldists[2,3],alldists[3,4],alldists[4,5])
  distnames = c("Lo.noneuc-Hi.noneuc","Hi.noneuc-Hi.noneuc","Hi.noneuc-Lo.noneuc","Lo.noneuc-Lo.noneuc")
  if(!is.null(covariate)) {
    dists=matrix(rep(dists,5),nrow=5)
    dists[2,] = c(lonu,hinu,hinu,lonu)
    dists[3,] = c(hinu,hinu,lonu,lonu)
    dists[4,] = c(locov,hicov,hicov,locov)
    dists[5,] = c(hicov,hicov,locov,locov)
    row.names(dists) = c("Distance","from.noneuc","to.noneuc",paste("from.",covariate,sep=""),paste("to.",covariate,sep=""))
    colnames(dists) = distnames
  } else {
    dists=matrix(rep(dists,3),nrow=3)
    dists[2,] = c(lonu,hinu,hinu,lonu)
    dists[3,] = c(hinu,hinu,lonu,lonu)
    row.names(dists) = c("Distance","from.noneuc","to.noneuc")
    colnames(dists) = distnames
  }
  dists
}

#' @title Least-cost path plot
#'
#' @description
#'  Plots 
#'  Requires packages secr, raster, gdistance, fields, igraph
#'  
#' @param mask A mask with covariate `noneuc` (obtained using \code{predictDsurface})
#' @param transitionFunction A user-defined function for calculating conductance between points 
#' (see function \code{transition} of the \code{gdistance} package).
#' 'userdfn1' in the \code{secr} vignette 'secr-noneuclidean.pdf', for example.)
#' @param covariate The name (a character) of a covariate to plot as an image plot
#' (typically because 'noneuc' depends on it).
#' @param n number of least-cost paths to plot
#' @param linecol color of least cost path line
#' @param directions parameter of function \code{transition}
#' @param symm parameter of function \code{transition}
#' @param directed parameter of function \code{transition}
#' @param ... arguments to function \code{plotcovariate}, which is used to plot the image.
#' 
#' @examples 
#' lcpathplot(Dhat.Nemegt,"stdBC",userdfun=function(x) 1/x,n=5,contour=FALSE,col=c("orange","brown"),key=FALSE, linecol="white",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1)

lcpathplot = function(mask,transitionFunction,covariate="noneuc",n=1,linecol="white", 
                      directions=16, symm=TRUE,directed=FALSE,...) {
  
  if(!is.element("noneuc",names(covariates(mask))))
    stop("Must have 'noneuc' as one of the mask covariates. It is not there.")
  if(!is.element(covariate,names(covariates(mask))))
    stop(paste("Invalid covariate: '",covariate,"'"," is not the name of one of the mask covariates.",sep=""))
  rastermask = raster(mask,"noneuc") # make raster with covariates(mask)$noneuc as values of pixels
  
  transfun=match.fun(transitionFunction)
  
  coords = coordinates(rastermask) # lookup table for vertex coordinates
  # secr models conductance (sigma) with a log link. Line below assumes that $noneuc is on linear predictor scale 
  # tr1<-transition(rastermask,transitionFunction=function(x) exp(lp(x)),directions=directions,symm=symm)
  tr1<-transition(rastermask,transitionFunction=transfun,directions=directions,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  #costs1<-costDistance(tr1CorrC,pts)
  
  plotcovariate(mask,covariate,...)
  
  dists = rep(NA,n) # to keep least-cost path distances in
  
  for(i in 1:n) {
    fromind = nearesttrap(unlist(locator(1)), mask)
    toind = nearesttrap(unlist(locator(1)), mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    to = c(x=mask$x[toind],y=mask$y[toind])
    from=matrix(from,ncol=2)
    to=matrix(to,ncol=2)
    npts=dim(from)[1]
    nptsto=dim(to)[1]
    if(nptsto != npts) stop("Must have same number of points in 'from' as in 'to'.")
    if(npts>1) pts = closest_coords(from[i,],to[i,],rastermask)
    else pts = closest_coords(from,to,rastermask)
    vpts = get_vertices(pts,rastermask)
    
    trmat=summary(transitionMatrix(tr1CorrC))
    #cbind(trmat,1/trmat$x)
    #    rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
    rel=data.frame(from=trmat$i,to=trmat$j,weight=trmat$x)
    #rel
    g = graph_from_data_frame(rel,directed=directed,vertices=NULL)
    #    attributes(g)$noneuc=1/trmat$x
    #    E(g)$weight=1/trmat$x
    attributes(g)$noneuc=trmat$x
    E(g)$weight=trmat$x
    #vertices = as_data_frame(g, what="vertices")
    #edges = as_data_frame(g, what="edges")
    svert=which(names(V(g))==vpts[1])
    evert=which(names(V(g))==vpts[2])
    spath=as.numeric(names(shortest_paths(g,from=svert,to=evert,weights=E(g)$weight)$vpath[[1]]))
    dists[i]=igraph:::distances(g,v=svert,to=evert,weights=attributes(g)$noneuc)
    
    nppts=length(spath)
    segments(coords[spath[-nppts],1],coords[spath[-nppts],2],coords[spath[-1],1],coords[spath[-1],2],col=linecol,lwd=lwd)
    points(coords[spath[c(1,nppts)],],pch=19,col="white",cex=1.5)
    points(coords[spath[c(1,nppts)],],pch=19,col=c("green","red"),cex=0.75)
  }
  
  invisible(dists)
}

#' @title Least-cost distance summary
#'
#' @description
#'  Tabulates least-cost distance between smallest and largest `noneuc` values on mask
#'  Requires packages secr, raster, gdistance, fields, igraph
#'  
#' @param masknu A mask with covariate `noneuc` (obtained using \code{predictDsurface})
#' @param userfn A user-defined function for calculating least-cost distances (see function 
#' 'userdfn1' in the \code{secr} vignette 'secr-noneuclidean.pdf', for example.)
#' @param covariate The name (a character) of a covariate to summarise with the distances 
#' (typically because 'noneuc' depends on it).
#' 
noneucdists = function(masknu,userfn,covariate=NULL) {
  if(!is.element("noneuc",names(covariates(masknu)))) stop("Need a mask covariate called `noneuc'.")
  if(!is.character(covariate)) stop("Argument `covariate' must be a character.")
  if(!is.null(covariate) & !is.element(covariate,names(covariates(masknu)))) 
    stop(paste("No covariate called",covariate,"in the mask"))
  f = match.fun(userfn)
  uNU = range(covariates(masknu)$noneuc)
  loNUind = which(covariates(masknu)$noneuc==uNU[1])[1]
  hiNUind = which(covariates(masknu)$noneuc==uNU[2])[1]
  lonu = covariates(masknu)$noneuc[loNUind]
  hinu = covariates(masknu)$noneuc[hiNUind]
  if(!is.null(covariate)) {
    covcol = which(names(covariates(masknu))==covariate)
    locov = covariates(masknu)[,covcol][loNUind]
    hicov = covariates(masknu)[,covcol][hiNUind]
  }
  # make mask with two points, distance 2 apart
  dets = trap.builder(frame=data.frame(x=c(-2,2),y=c(0,0)),method="all")
  msk = make.mask(dets,buffer=2, spacing=1)
  # calculate dists for lo-lo, lo-hi, hi-lo, hi-hi
  covariates(msk) = data.frame(noneuc=rep(c(lonu,hinu,hinu,lonu,lonu,rep(0,dim(msk)[1]-5))))
  alldists = f(as.matrix(msk),as.matrix(msk),msk)
  dists = c(alldists[1,2],alldists[2,3],alldists[3,4],alldists[4,5])
  distnames = c("Lo.noneuc-Hi.noneuc","Hi.noneuc-Hi.noneuc","Hi.noneuc-Lo.noneuc","Lo.noneuc-Lo.noneuc")
  if(!is.null(covariate)) {
    dists=matrix(rep(dists,5),nrow=5)
    dists[2,] = c(lonu,hinu,hinu,lonu)
    dists[3,] = c(hinu,hinu,lonu,lonu)
    dists[4,] = c(locov,hicov,hicov,locov)
    dists[5,] = c(hicov,hicov,locov,locov)
    row.names(dists) = c("Distance","from.noneuc","to.noneuc",paste("from.",covariate,sep=""),paste("to.",covariate,sep=""))
    colnames(dists) = distnames
  } else {
    dists=matrix(rep(dists,3),nrow=3)
    dists[2,] = c(lonu,hinu,hinu,lonu)
    dists[3,] = c(hinu,hinu,lonu,lonu)
    row.names(dists) = c("Distance","from.noneuc","to.noneuc")
    colnames(dists) = distnames
  }
  dists
}


#' @title Plots least cost paths on mask
#'
#' @description
#'  Plots least cost distance between pairs of points, together (optionally) with mask cost surface. Marks start with 
#'  a green dot, and end with a red dot, and joins them by line segments comprising the least cost path.
#'  
#'  Requires function calc.p.logPis() and plotcovariate().
#'  Requires packages secr, raster, gdistance, fields, igraph
#'  
#' @param from matrix or data frame of coordinates from which to start, with an (x,y) pair per line
#' @param to matrix or data frame of coordinates at which to end, with an (x,y) pair per line 
#' dim(from) must be equal to dim(to)
#' @param mask \code{secr} mask object. Must have covariate called 'noneuc' containing cost
#' @param costname Name of variable to use in cost calculation
#' @param costfun Cost function name
#' @param directed If TRUE, use directed graph for transition between cells, else undirected 
#' @param symm If TRUE, cost is same in both directions 
#' @param plotmask If TRUE, plots mask cost surface
#' @param add  If TRUE, adds to an existing plot
#' @param lincol Line color for least cost path
#' @param lwd line width for least cost path
#' 
#' @examples 
#' library(gdistance)
#' library(igraph)
#' library(fields)
#' library(raster)
#' library(secr)
#' 
#' ny=4; nx=4 # dimensions of mask
#' # set costs (NA is "hole" - nothing there & can't go there):
#' costs=c(100,100,100,100,NA,100,100,NA,1,NA,100,1,1,1,1,1) 
#' rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),noneuc=costs) # make data frame with coords and costs
#' 
#' rmask=read.mask(data=rmesh,columns="noneuc")  # make mask with covariate 'noneuc' containing cost
#' plotcovariate(rmask,covariate="noneuc",contour=FALSE) # look at the covariate surface
#' 
#' from=c(4.1,0.9); to=c(2.3,3.8) # corrdinates of start and end points (Note: need not be exactly on mask coords)
#' from=data.frame(x=c(4.1,1,2),y=c(0.9,1,2)); to=data.frame(x=c(2,2,2),y=c(4,4,4))
#' dist = plot_lcpath(from,to,rmask,lwd=2,linecol="white") # calculate and plot the least cost path, returning cost
#' dist # print the cost
#' 
#' cfun=function(x) exp(diff(x)) # asymmetric cost function
#' 
#' plot_lcpath(from,to,rmask,costfun=cfun,lwd=2,linecol="white") # use different cost function
#' dist # print the cost
#' 
#' cfun=function(x) exp(diff(x)) # asymmetric cost function
#' 
#' plot_lcpath(from,to,rmask,costfun=cfun,lwd=2,linecol="white") # use different cost function
#' 
#' plot_lcpath(from,to,rmask,costfun=cfun,symm=FALSE,lwd=2,linecol="white") # without symmetry, undirected
#' 
#' plot_lcpath(from,to,rmask,costfun=cfun,symm=FALSE,directed=TRUE,lwd=2,linecol="white") # without symmetry, directed
#' 
#' plot_lcpath(from,to,rmask,costfun=cfun,lwd=2,linecol="white") # with symmetry, undirected
#' 
#' @export plot_lcpath
#' 
plot_lcpath = function(from,to,mask,costname="noneuc",costfun="mean",plotcovariate="noneuc",directed=FALSE,symm=FALSE,
                       plotmask=TRUE,add=TRUE,linecol="black",lwd=1,col=tim.colors(40),key=TRUE) {
  
  
  if(!is.element(costname,names(covariates(mask))))
    stop(paste("Invalid costname: '",costname,"'"," is not the name of one of the mask covariates.",sep=""))
  if(!is.element(plotcovariate,names(covariates(mask))))
    stop(paste("Invalid plotcovariate: '",plotcovariate,"'"," is not the name of one of the mask covariates.",sep=""))
  rastermask = raster(mask,costname) # make raster with covariates(mask)$costname as values of pixels
  
  f=match.fun(costfun)
  
  coords = coordinates(rastermask) # lookup table for vertex coordinates
  tr1<-transition(rastermask,transitionFunction=function(x) 1/f(x),directions=8,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  #costs1<-costDistance(tr1CorrC,pts)
  
  if(is.vector(from)) from=matrix(from,ncol=2)
  if(is.vector(to)) to=matrix(to,ncol=2)
  npts=dim(from)[1]
  nptsto=dim(to)[1]
  if(nptsto != npts) stop("Must have same number of points in 'from' as in 'to'.")
  for(i in 1:npts) {
    if(npts>1) pts = closest_coords(from[i,],to[i,],rastermask)
    else pts = closest_coords(from,to,rastermask)
    vpts = get_vertices(pts,rastermask);vpts
    
    trmat=summary(transitionMatrix(tr1CorrC))
    #cbind(trmat,1/trmat$x)
    rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
    #rel
    g = graph_from_data_frame(rel,directed=directed,vertices=NULL)
    attributes(g)$noneuc=1/trmat$x
    E(g)$weight=1/trmat$x
    #vertices = as_data_frame(g, what="vertices")
    #edges = as_data_frame(g, what="edges")
    svert=which(names(V(g))==vpts[1])
    evert=which(names(V(g))==vpts[2])
    spath=as.numeric(names(shortest_paths(g,from=svert,to=evert,weights=E(g)$weight)$vpath[[1]]))
    dist=igraph:::distances(g,v=svert,to=evert,weights=attributes(g)$noneuc)
    
    nppts=length(spath)
    #as.matrix(rastermask)
    if(i==1) {
      if (plotmask) plotcovariate(mask,covariate=plotcovariate,contour=FALSE,col=col,key=key,asp=1)
      if(!add) plot(coords,type="n")
    }
    segments(coords[spath[-nppts],1],coords[spath[-nppts],2],coords[spath[-1],1],coords[spath[-1],2],col=linecol,lwd=lwd)
    points(coords[spath[c(1,nppts)],],pch=19,col="white",cex=1.5)
    points(coords[spath[c(1,nppts)],],pch=19,col=c("green","red"),cex=0.75)
  }
  
  
  invisible(dist)
  
}


#' @title Finds closest coordinates on raster to two points
#'
#' @description
#'  Uses function over() from package sp to overlay points on raster and return closest raster coordinates
#'  
#' @param from pair of coordinates (x,y) from which to start
#' @param to pair of coordinates (x,y) to get to
#' @param rastermask Raster object (typically created from mask by something like 
#' rastermask = raster(mask,"noneuc"))
#' 
#' @return Returns the coordinates of teh closest point on the raster, as a matrix with two columns (x,y), 
#' named s1 and s2.
closest_coords=function(from,to,rastermask){
  ends=SpatialPoints(rbind(from,to))
  grid=as(rastermask, 'SpatialGrid') 
  xy=over(ends,grid)
  return(coordinates(grid)[xy,])
}


#' @title Finds vertex index on graph made from raster
#'
#' @description
#'  Finds vertex index on graph made from raster, of points at coordinates pts. Vertex index is just the row of 
#'  the point in the raster object.
#'  
#' @param pts Matrix whose rows are (x,y) coordinates of points on raster
#' @param raster Raster object.
#' 
#' @return Returns the row numbers of raster that correpond to pts. Note that pts must match exactly some 
#' coordinates of raster (use \code{closest_coords} to find closest coordinates if necessary).
#' 
get_vertices = function(pts,rastermask){
  coords = coordinates(rastermask) # lookup table from index produced by transition() to coordinates
  npts = dim(pts)[1]
  vert = rep(NA,npts)
  for(i in 1:npts){
    vert[i] = which(coords[,1]==pts[i,1] & coords[,2]==pts[i,2])
  }
  return(vert)
}

#' @title Creates the igraph of a mask object
#'
#' @description
#'  Creates an igraph object with a vertex for each mask point and edges to neighbours, weighted according 
#'  to the cost function \code{costfun}, using the mask covariate \code{costname}.
#'  
#'  Requires packages raster, gdistance, igraph
#'  
#' @param mask \code{secr} mask object. Must have covariate called 'noneuc' containing cost
#' @param costname Name of variable to use in cost calculation
#' @param costfun Cost function name
#' @param directed If TRUE, use directed graph for transition between cells, else undirected 
#' @param symm If TRUE, cost is same in both directions 
#' 
#' @examples 
#' ny=4; nx=4 # dimensions of mask
#' # set costs (NA is "hole" - nothing there & can't go there):
#' costs=c(100,100,100,100,NA,100,100,NA,1,NA,100,1,1,1,1,1) 
#' rmesh=data.frame(x=rep(1:nx,ny),y=rep(1:ny,rep(nx,ny)),noneuc=costs) # make data frame with coords and costs
#' 
#' rmask=read.mask(data=rmesh,columns="noneuc")  # make mask with covariate 'noneuc' containing cost
#' ig=make_igraph(rmask,"noneuc")
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' cfun=function(x) exp(diff(x)) # asymmetric cost function
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun,directed=TRUE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun,directed=TRUE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
#' ig=make_igraph(rmask,"noneuc",costfun=cfun,directed=TRUE,symm=FALSE)
#' plot(ig, edge.label=round(E(g)$weight, 3))
#' 
make_igraph = function(mask,costname,costfun="mean",directed=FALSE,symm=TRUE) {
  
  require(raster)
  require(gdistance)
  
  if(!is.element(costname,names(covariates(mask))))
    stop(paste("'",costname,"'"," is not the name of one of the mask covariates.",sep=""))
  rastermask = raster(mask,costname) # make raster with covariates(mask)$costname as values of pixels
  
  f=match.fun(costfun)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/mean(x),directions=8)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/exp(diff(x)),directions=8,symm=FALSE)
  tr1<-transition(rastermask,transitionFunction=function(x) 1/f(x),directions=8,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  #costs1<-costDistance(tr1CorrC,pts)
  
  pts = closest_coords(from,to,rastermask);pts
  vpts = get_vertices(pts,rastermask);vpts
  
  trmat=summary(tr1CorrC)
  rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
  if(directed) g = graph_from_data_frame(rel,directed=TRUE,vertices=NULL)
  else g = graph_from_data_frame(rel,directed=FALSE,vertices=NULL)
  attributes(g)$noneuc=1/trmat$x
  E(g)$weight=1/trmat$x
  
  return(g)  
}
