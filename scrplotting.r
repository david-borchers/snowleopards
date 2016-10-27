#' @title Plots overall detection probability surface p..(s) on mesh
#'
#' @description
#'  Plots overall detection probability surface p..(s) on mesh. Currenly only plots for a 
#'  single individual and only multi-catch traps and count detectors with a single session 
#'  implemented.
#'  
#'  Requires function calc.p.logPis() and plotcovariate().
#'  
#' @param capthist Object of class  '\code{capthist}' from package \code{secr}. Currently only 
#' multi-catch traps and count detectors with a single session implemented.
#' @param mesh Object of class  '\code{capthist}' from package \code{secr}.
#' @param pars Vector of parameters, in the format required by function \code{scr.negloglik}.
#' @param dist K by M Matrix distances from each of the K detectors to each of the M points on 
#' the mesh.
#' @param contour TRUE if you want a contour plotted on top of the image plot.
#' @param key TRUE if you want a key (legend) added showing probability levels.
#' @param ... Other arguments for \code{plotcovariate}.
#' 
#' @export plot.p..
plot.p..=function(capthist, mesh, pars, dist=NULL, contour=TRUE, key=TRUE, ...) {
  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.logPis(capthist, mesh, pars, dist=dist)
  
  covariates(mesh)$p..s = p..s
  plotcovariate(mesh,covariate="p..s",main="p..s",contour=contour, key=key, ...)
}



#' @title Image plot contour of activity centre location.
#'
#' @description
#'  Image plot contour of activity centre location. Currenly only plots for a single individual
#'  and only multi-catch traps and count detectors with a single session implemented.
#'  
#'  Requires function calc.p.logPis() and plotcovariate().
#'  
#' @param i A scalar or vector indices of individuals to plot - corresponding to the order they 
#' appear in the capture history.
#' @param capthist Object of class  '\code{capthist}' from package \code{secr}. Currently only 
#' multi-catch traps and count detectors with a single session implemented.
#' @param mesh Object of class  '\code{capthist}' from package \code{secr}.
#' @param pars Vector of parameters, in the format required by function \code{scr.negloglik}.
#' @param dist K by M Matrix distances from each of the K detectors to each of the M points on 
#' the mesh.
#' @param contour TRUE if you want a contour plotted on top of the image plot.
#' @param key TRUE if you want a key (legend) added showing probability levels.
#' @param ... Other arguments for \code{plotcovariate}.
#' 
#' @export plot.Pi
plot.Pi=function(i,capthist, mesh, pars, dist=NULL, contour=TRUE, key=TRUE, ...) {
  M=dim(mesh)[1] # number of mesh points
  a = rep(attributes(mesh)$area,M) # area of each mesh cell (vector of length M)
  
  dets=traps(capthist) # detectors
  
  D = rep(exp(pars["D"]),M) # Density at each mesh point (vector of length M)
  
  # Calculate p..(s) and log(P_i(s)) at each mesh point
  scr = calc.p.logPis(capthist, mesh, pars, dist=dist)
  
  covariates(mesh)$Pi.s = scr$Pi.s[i,]/sum(scr$Pi.s[i,])
  plotcovariate(mesh,covariate="Pi.s",main=expression(f(s[i]*"|"*capthist[i])),contour=contour, key=key, ...)
  plot(dets,add=TRUE,detpar=list(pch=19,col="black"))
  if(dim(ch)[2]==1) freq = apply(apply(ch,c(2,3),sum),2,sum)
  else freq=apply(capthist[i,,],2,sum)
  detected=which(freq>0)
  text(dets$x[detected],dets$y[detected],labels=freq[detected],cex=0.75,col="white")
}


#' @title Draws histogram.
#'
#' @description
#'  Utility function to draw histograms with more options than \code{hist} allows.
#'  
#' @param height Height of histogram bars.
#' @param breaks Locations of boundaries of histogram bins (must be 1 longer than \code{height}).
#' @param lineonly If TRUE, uses \code{\link{lines}} to draw lines on current plot; else uses 
#' \code{\link{plot}} to draw lines on new plot.
#' @param outline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param fill If TRUE, uses polygon() to fill barsl in this case valid arguments to polygon() 
#' are passed via argument(s) "...". If fill==FALSE, valid arguments to plot() or lines() are 
#' passed via argument(s) "..."
#' @param ylim Range of y-axis.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param ... See aargument \code{fill}.
histline=function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,ylim=range(height),
                  xlab="x",ylab="y",...)
{
  n=length(height)
  if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
  if(outline) {
    y=c(0,rep(height,times=rep(2,n)),0)
    x=rep(breaks,times=rep(2,(n+1)))
  }   else {
    y=rep(0,4*n)
    x=rep(0,4*n+2)
    for(i in 1:n) {
      y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
      x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
    }
    x=x[1:(4*n)]
  }
  if(lineonly) {
    if(!fill) lines(x,y,...)
    else polygon(x,y,...)
  } else {
    if(!fill) plot(x,y,type="l",ylim=ylim,xlab=xlab,ylab=ylab,...)
    else {
      plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab)
      polygon(x,y,...)
    }
  }
}



#' @title Encounter rate fit plot
#'
#' @description
#'  Plots mean encounter rate by distance interval, with estimated encounter rate function
#'  overlaid. NB: Distances calculated from estimated activity centre mode. Don't over-interpret
#'  the plot - the activity centre locations may have very high uncertainty attached to them, 
#'  and hence the distances may be very uncertain!
#'  
#'  The function only works with traps of class '\code{count}'.
#'  
#'  Uses function '\code{histline}'.
#'  
#' @param scrfit Object of class  '\code{secr}' from package \code{secr}. The function only 
#' works with traps of class '\code{count}'.
#' @param dmax Maximum distance to plot.
#' @param binwidth width of each distance interval to plot.
#' @param lwd Line width for encounter rate function.
#' @param xlab Label for x-axis.
#' @param ylab Label for y-axis.
#' @param hcol Colour of histogram.
#' @param hlineonly If TRUE, uses \code{\link{lines}} to draw histogram lines on current plot; 
#' else uses \code{\link{plot}} to draw lines on new plot.
#' @param houtline If TRUE, draws only the outline (profile) of the histogram; else draws each 
#' complete bar.
#' @param hfill TRUE, fills bards with \code{hcol}, else leaves them empty.
#' 
#' @export plot.er.fit
plot.er.fit=function(scrfit,dmax,binwidth,lwd=2,xlab="x",ylab="y",
                     hcol="red",hlineonly=FALSE,houtline=FALSE,hfill=TRUE) {
  
  if (detector(traps(scrfit$capthist))!="count") stop("This function only works with count detectors.")
  
  ch = scrfit$capthist
  n=summary(ch)$counts["n","Total"] # number of individuals
  loc=data.frame(x=rep(0,n),y=rep(0,n)) # data frame for locations of individuals
  for(i in 1:n) loc[i,] = fxi.mode(er.fit,i=i) # add estimated modes of invivs to data frame
  dists = distances(loc,kruger.cams) # distances between indiv modes and detectors
  capts = apply(kruger.capt,c(1,3),sum) # sum over occasions
  dvec = as.vector(dists) # turn into a vector of distances
  cvec = as.vector(capts) # turn into a vector of counts
  keep = dvec<=dmax # look only within given max dist
  
  breaks=seq(0,dmax,binwidth)
  nb=length(breaks)-1
  mean.n=rep(0,nb)
  for(i in 1:nb) {
    keep=which(breaks[i]<=dvec & dvec<breaks[i+1])
    mean.n[i]=mean(cvec[keep])
  }
  
  # gymnastics to avoid plotting while getting plot values in pl
  ff <- tempfile()
  png(filename=ff)
  pl=plot(scrfit,xval=seq(0,dmax,length=200),lwd=2)
  dev.off()
  unlink(ff)
  
  # now do the plot we want
  histline(mean.n,breaks,xlab=xlab,ylab=ylab,col=hcol,fill=hfill,lineonly=hlineonly,outline=houtline,
           ylim=c(0,max(c(mean.n,pl$y))))
  plot(scrfit,xval=seq(0,dmax,length=200),lwd=lwd,add=TRUE)
}



#' @title Plots image (and optionally contours) of mask covariate value
#' @param mask is an object of class `mask'
#' @param covariate is a character variable with the name of one of the covariates in mask (the one to plot)
#' @param contour is a logical, TRUE if want contour plots on image
#' @param ... other arguments to be passed to \code{prep4image}
#' @export
plotcovariate=function(mask,covariate,contour=TRUE,key=TRUE, ...) {
  cnum=which(names(covariates(mask))==covariate)
  if(is.null(cnum)) stop(paste("No covariate(s) called",covariate))
  if(length(cnum)>1) warning("Can only plot one covariate at a time. First covariate being plotted.")
  dat=data.frame(x=mask$x,y=mask$y,z=covariates(mask)[[cnum]])
  prep4image(dat,contour=contour,key=key,...)
}

#' @title Plots image (and optionally contours) of density from Dsurface object
#' @param Dsurface is an object of class `Dsurface'
#' @param contour is a logical, TRUE if want contour plots on image
#' @param ... other arguments to be passed to \code{prep4image}
#' @export
plot.Dsurface=function(Dsurface,contour=TRUE,key=TRUE, ...) {
  dat=data.frame(x=Dsurface$x,y=Dsurface$y,z=covariates(Dsurface)$D.0)
  prep4image(dat,contour=contour,key=key,...)
}


#' @title Prepares data frame for plotting with image/contour/persp.
#'   
#' @description From an input data frame with columns x, y and z, this function 
#'   creates a list with elements x, y and z in a format suitable for passing to
#'   functions \code{\link{image}}, \code{\link{contour}} or 
#'   \code{\link{persp}}. The coordinates in \code{data} are assumed to come
#'   from a 2D grid of points.
#'   
#' @param data a data frame with columns x and y being Cartesian coordinates, 
#'   and z being the values of some variable at each coordinate.
#' @param plot if \code{TRUE} then an image plot will be drawn using 
#'   \code{\link{image.plot}}
#' @param contour if \code{TRUE} then contours will be added (only used when 
#'   \code{plot=TRUE})
#' @param key logical for whether or not to include key when \code{plot = TRUE} (\code{\link{image.plot}} is used when \code{key = TRUE}, \code{\link{image}} is used when \code{key = FALSE})
#' @param ... other arguments to pass to \code{\link{image}} or \code{\link{image.plot}} (only used 
#'   when \code{plot=TRUE})
#'   
#' @details Sorts z on values of x first, then y, then creates a matrix of 
#'   z-values from this. Returns a list with elements x (unique values of x, in 
#'   increasing order), y (unique values of y, in increasing order) and z 
#'   (matrix of z-values in appropriate order for image/contour/persp). 
#'   
#'   If the original z is a factor variabele, the z returned is a matrix of integers 
#'   between 1 and length(levels(z)) and the output list has an attributes called 
#'   ``facnames'' that is a character vector containing the levels as factor 
#'   variables, with z=1 corresponding to the first name, z=2 to the second, etc.
#' @export
prep4image = function(data, plot = TRUE, contour = TRUE, key = TRUE, ...){
  
  # convert factor data$z to integer:
  zfactor=FALSE
  if(is.factor(data$z)) {
    zfactor=TRUE
    fac=data$z
    facnames=levels(fac)
    nlevels=length(facnames)
    data$z=rep(0,length(fac))
    got=rep(FALSE,nlevels)
    for(i in 1:nlevels){
      j=which(fac==facnames[i])
      if(length(j)>0) got[i]=TRUE
      data$z[j]=i
    }
    facnames=facnames[got] # remove factor names not in mask
  }
  data = as.matrix(data)
  
  x = sort(unique(data[,"x"]))
  y = sort(unique(data[,"y"]))
  
  z = matrix(NA, nrow = length(x), ncol = length(y))
  
  for(i in 1:length(x)){
    for(j in 1:length(y)){
      m = which(data[,"x"] == x[i] & data[,"y"] == y[j]) ; m
      z[i,j] = if(length(m) == 0) NA else data[,"z"][m]
    }
  }
  
  if(plot){
    if(key){
      image.plot(x, y, z, ...)
    }else{
      image(x, y, z, ...)
    }
    if(contour) contour(x, y, z, add = TRUE)
  }
  
  outlist=list(x = x, y = y, z = z)
  if(zfactor) attributes(outlist)$facnames=facnames
  
  invisible(outlist)
  
}


#' @title add new covariate consisting of distances to points that have covariate==distance.to
#' @param mask is object of class "mask".
#' @param covariate is name of one of the elements of covariates(mask)
#' @param distance.to is the target value of covariate; distances to mask points with this value are calculated.
#' @param dname is the name you want to give to the new covariate
add.dist.to=function(mask,covariate,distance.to,dname=NULL,overwrite=FALSE){
  covs=names(covariates(mask))
  if(is.null(covs)) stop("No covariates in mask. Can't add distance to anything")
  ncov=length(covs)
  if(is.null(dname)) dname=paste("dist2",covariate,"=",distance.to,sep="")
  if(is.element(covariate,covs)) {
    if(dname %in% names(covariates(mask))) {
      if(overwrite) {
        warning(paste("Covariate called",dname,"already existed. It has been overwritten."))
        covi=which(names(covariates(mask))==dname)
      } else stop(paste("Covariate called",covariate,"already exists. Use overwrite=TRUE if you want to overwrite it."))
    } else {
      covi=ncov+1
    }
    cov=covariates(mask)[which(covs==covariate)]
    if(is.null(cov)) stop(paste("No covariate called",covariate,"in mask."))
    targets=which(cov==distance.to)
    distances=rdist(mask,mask[targets,])
    d2=apply(distances,1,min)
    covariates(mask)[[covi]]=d2
    names(covariates(mask))[[covi]]=dname
    return(mask)
  } else stop(paste("No covariate called",covariate,"in mask."))
}
