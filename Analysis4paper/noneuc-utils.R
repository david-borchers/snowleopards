geommean = function(x) exp(mean(x))

geomLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = geommean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

arithmean = function(x) mean(exp(x))

arithLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = arithmean, directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}


lcusageplot = function(fit,n=512,mask=NULL,base="noneuc.0",lcdfun="geomLCdist",...) {
  if(is.null(mask)) mask = fit$mask
  
  lcd.fun = match.fun(lcdfun)
  
  if(!is.element("noneuc.0",names(covariates(mask))))
    stop("Must have 'noneuc.0' as one of the mask covariates. It is not there.")
  if(!is.element(base,names(covariates(mask)))) {
    warning(paste("mask does not have a covariate called ",base,"; noneuc.0 being used instead."))
  }
  covariates(mask)$base = covariates(mask)$noneuc.0
  
  for(i in 1:n){
    plotcovariate(mask,"base",col=parula(40),what="image")
    fromind = nearesttrap(unlist(locator(1)), mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    from=matrix(from,ncol=2)
    dists = lcd.fun(from,mask,mask)
    dfn <- secr:::getdfn(fit$detectfn)
    p = dfn(dists[1,],unlist(detectpar(fit)))
    covariates(mask)$p = p/sum(p)
    plotcovariate(mask,"p",what="image",...)
    points(from,col="white",pch=19)
    points(from,cex=0.8,pch=19)
    waitclick = unlist(locator(1))
  }
  
  invisible(mask)
}

lcpathplot = function(mask,transitionFunction,type="noneuc",n=512,background="d2.river",linecol="white", 
                      directions=16, symm=TRUE,directed=FALSE,lwd=1,...) {
  
  if(!is.element("noneuc.0",names(covariates(mask))))
    stop("Must have 'noneuc.0' as one of the mask covariates. It is not there.")
  if(!is.element(type,c("noneuc","sigma")))
    stop(paste("Invalid type: '",type,"'"," It must be `noneuc` or `sigma`.",sep=""))
  rastermask = raster(mask,"noneuc.0") # make raster with covariates(mask)$noneuc.0 as values of pixels
  
  transfun=match.fun(transitionFunction)
  
  coords = coordinates(rastermask) # lookup table for vertex coordinates
  # secr models conductance (sigma) with a log link. Line below assumes that $noneuc.0 is on linear predictor scale 
  # tr1<-transition(rastermask,transitionFunction=function(x) exp(lp(x)),directions=directions,symm=symm)
  tr1<-transition(rastermask,transitionFunction=transfun,directions=directions,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  
  if(type=="noneuc") plotcovariate(mask,"noneuc.0",...)
  else {
    covariates(mask)$sigma = exp(covariates(mask)$noneuc.0)
    plotcovariate(mask,"sigma",...)
  }
  
  dists = rep(NA,n) # to keep least-cost path distances in
  
  for(i in 1:n) {
    fromind = nearesttrap(unlist(locator(1)), mask)
    toind = nearesttrap(unlist(locator(1)), mask)
    from = c(x=mask$x[fromind],y=mask$y[fromind])
    to = c(x=mask$x[toind],y=mask$y[toind])
    from=matrix(from,ncol=2)
    to=matrix(to,ncol=2)
    #    npts=dim(from)[1]
    #    nptsto=dim(to)[1]
    #    if(nptsto != npts) stop("Must have same number of points in 'from' as in 'to'.")
    #    if(npts>1) pts = closest_coords(from[i,],to[i,],rastermask)
    #    else pts = closest_coords(from,to,rastermask)
    pts = closest_coords(from,to,rastermask)
    vpts = get_vertices(pts,rastermask)
    
    trmat=summary(transitionMatrix(tr1CorrC))
    #cbind(trmat,1/trmat$x)
    #    rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
    rel=data.frame(from=trmat$i,to=trmat$j,weight=trmat$x)
    #rel
    g = graph_from_data_frame(rel,directed=directed,vertices=NULL)
    #    attributes(g)$noneuc.0=1/trmat$x
    #    E(g)$weight=1/trmat$x
    attributes(g)$noneuc.0=trmat$x
    E(g)$weight=trmat$x
    #vertices = as_data_frame(g, what="vertices")
    #edges = as_data_frame(g, what="edges")
    svert=which(names(V(g))==vpts[1])
    evert=which(names(V(g))==vpts[2])
    # NB: Need to invert E(g) so that higher values lead to shorter distances:
    spath=as.numeric(names(shortest_paths(g,from=svert,to=evert,weights=1/E(g)$weight)$vpath[[1]]))
    dists[i]=igraph:::distances(g,v=svert,to=evert,weights=attributes(g)$noneuc.0)
    
    nppts=length(spath)
    segments(coords[spath[-nppts],1],coords[spath[-nppts],2],coords[spath[-1],1],coords[spath[-1],2],col=linecol,lwd=lwd)
    points(coords[spath[c(1,nppts)],],pch=19,col="white",cex=1.5)
    points(coords[spath[c(1,nppts)],],pch=19,col=c("green","red"),cex=0.75)
  }
  
  invisible(dists)
}

closest_coords=function(from,to,rastermask){
  ends=SpatialPoints(rbind(from,to))
  grid=as(rastermask, 'SpatialGrid') 
  xy=over(ends,grid)
  return(coordinates(grid)[xy,])
}




plotcovariate = function(mask,covariate, ...) {
  require(sp)
  
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  pts = cbind(mask$x,mask$y)
  spdfdat = data.frame(covariate=covariates(mask)[[cnum]])
  spdf = SpatialPixelsDataFrame(pts, spdfdat)
  plot(spdf, ...)
  
  invisible(spdf)
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
#' @return Returns the coordinates of the closest point on the raster, as a matrix with two columns (x,y), 
#' named s1 and s2, with first row corresponding to 'from' coordinates, and second row corresponding to 'to' 
#' coordinates.
#' @export closest_coords
#' 
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
#' @export get_vertices
#' 
get_vertices = function(pts,rastermask){
  target = nearest.raster.point(pts[,1],pts[,2],as.im(rastermask),indices=FALSE)
  coords = coordinates(rastermask) # lookup table from index produced by transition() to coordinates
  npts = dim(pts)[1]
  vert = rep(NA,npts)
  for(i in 1:npts){
    #    vert[i] = which(coords[,1]==target$x[i] & coords[,2]==target$y[i])
    dst = sqrt((coords[,1]-target$x[i])^2 + (coords[,2]-target$y[i])^2)
    vert[i] = which(dst == min(dst))[1]
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
#' @export make_igraph
#' 
make_igraph = function(mask,costname,costfun="mean",directed=FALSE,symm=TRUE) {
  
  
  if(!is.element(costname,names(covariates(mask))))
    stop(paste("'",costname,"'"," is not the name of one of the mask covariates.",sep=""))
  rastermask = raster(mask,costname) # make raster with covariates(mask)$costname as values of pixels
  
  f=match.fun(costfun)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/mean(x),directions=8)
  #tr1<-transition(rastermask,transitionFunction=function(x) 1/exp(diff(x)),directions=8,symm=FALSE)
  tr1<-transition(rastermask,transitionFunction=function(x) 1/f(x),directions=8,symm=symm)
  tr1CorrC<-geoCorrection(tr1,type="c",multpl=FALSE,scl=FALSE)
  #costs1<-costDistance(tr1CorrC,pts)
  
  pts = closest_coords(from,to,rastermask)
  vpts = get_vertices(pts,rastermask)
  
  trmat=summary(tr1CorrC)
  rel=data.frame(from=trmat$i,to=trmat$j,weight=1/trmat$x)
  if(directed) g = graph_from_data_frame(rel,directed=TRUE,vertices=NULL)
  else g = graph_from_data_frame(rel,directed=FALSE,vertices=NULL)
  attributes(g)$noneuc=1/trmat$x
  E(g)$weight=1/trmat$x
  
  return(g)  
}


plotcovariate = function(mask,covariate, ...) {
  require(sp)
  
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  pts = cbind(mask$x,mask$y)
  spdfdat = data.frame(covariate=covariates(mask)[[cnum]])
  spdf = SpatialPixelsDataFrame(pts, spdfdat)
  plot(spdf, ...)
  
  invisible(spdf)
}
