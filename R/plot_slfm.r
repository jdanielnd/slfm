#' Plot the data matrix along with probabilities box-plots
#'
#' This function is used to fit a Bayesian sparse
#' latent factor model to a directory of numeric matrices.
#'
#' @param y data matrix
#' @param result matrix of p_value chains from SLFM
#' @param stand indicates whether the data should be standardized
#' @param reordRow indicates whether the rows should be ordered
#' @param reordCol indicates whether the cols should be ordered
#' @param highContrast indicates if high contrast should be used
#' @importFrom lattice levelplot
#' @importFrom lattice bwplot
#' @importFrom reshape melt
## ' @export
plot_slfm <- function(
  y, result, clas, stand = TRUE, reordRow = TRUE,
  reordCol = TRUE, highContrast = TRUE) {

  y = as.matrix(y)
  xl = "matricesroarrays" 
  yl = "Probes"

  if(stand) {
    y = t(apply(y,1,scale)) # standardize the rows of y
  }

  if(reordRow) {
    medrow = apply(y,1,median) # median of the rows of y
    y = y[order(medrow),] # order the rows of y 
  }

  if(reordCol) {
    medcol = apply(y,2,median) # median of the columns of y
    y = y[,order(medcol)] # order the columns of y
  }

  if(highContrast) {
    y = (abs(y)^(1/3))*sign(y) # higher contrast
  }

  mi = as.numeric(quantile(y,0.001))
  ma = as.numeric(quantile(y,0.999))
  nr = nrow(y)
  nc = ncol(y)

  if (nr==nc) { # Image of Correlation and Covariance matrices
    xl=""
    yl=""
    if (ma==1) {
      mi = -1
    } else {
      mi=-ma
    }
  }
  if (nr<=5) {
    spr=1
  } else {
    spr=round(nr/10)
  }
  if (nc<=5) {
    spc=1
  } else {
    spc=round(nc/10)
  }
  sc = list(
    x=list(
      at=c(seq(1,nc, spc),nc),
      labels=c(seq(1,nc,spc),nc)
    ), 
    y=list(
      draw=FALSE
    )
  )
  col.l = colorRampPalette(c('blue','white','red'))
  cbar = seq(mi,ma,length.out=100)

  ord = order(apply(result, 2, median))

  theme.novpadding <-
     list(layout.widths =
          list(left.padding = 0,
        key.ylab.padding = 0,
        ylab.axis.padding = 0
  ))

  lp = levelplot(t(y[ord,]),col.regions=col.l,xlab=xl,ylab="",scales=sc, 
    at=cbar,aspect="fill", par.settings = theme.novpadding)

  colnames(result) <- formatC(1:ncol(result), width=round(log(ncol(result),10)+1), flag="0")
  res_melt = reshape::melt(result)[,-1]
  colnames(res_melt) = c("mat","val")
  res_melt[,1] <- factor(res_melt[,1], levels=ord)

  scbp = list( x=list(at=c(0,0.5,1),labels=c(0,0.5,1)), 
          y=list(at=1:nr,labels=(1:nr)[ord]) )

  theme.novpadding2 <-
    list(
      layout.widths = list(
        left.padding = 0,
        key.ylab.padding = 0,
        ylab.axis.padding = 0,
        right.padding = 0
      )
    )

  bp = bwplot(mat~val, data=res_melt, horizontal = TRUE,
    as.table = TRUE, xlab = "p(z_i=1|-)", scales = scbp,
    par.settings = theme.novpadding2)

  theme.novpadding3 <-
   list(layout.widths =
        list(left.padding = 0,
      key.ylab.padding = 0,
      ylab.axis.padding = 0,
      axis.key.padding = 0,
      right.padding = 0))

  sccp = list( x=list(at=NULL,labels=NULL), 
          y=list(at=1:nr,labels=(1:nr)[order(medrow)]) )
  col.l2 = colorRampPalette(c('white','black'))

  cp = levelplot(t(as.numeric(clas[order(medrow)])), colorkey = NULL, par.settings=theme.novpadding3,
    scales=sccp, xlab = "Call", col.regions=col.l2)
  
  print(cp,position=c(0,0,0.1,1),more=T)  
  print(bp,position=c(0.1,0,0.25,1), more=T)
  print(lp,position=c(0.25,0,1,1))

}