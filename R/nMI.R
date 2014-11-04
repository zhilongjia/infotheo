#' normalised mutual information computation
#' 
#' \code{nMI} takes two random variables as input and computes the mutual 
#' information in nats according to the entropy estimator \code{method}.If Y is 
#' not supplied and X is a matrix-like argument, the function returns a matrix 
#' of mutual information between all pairs of variables in the dataset X. 
#' 
#' method of entropy estimator:
#' \itemize{
#' \item "emp" : This estimator computes the entropy of the empirical probability distribution.
#' \item "mm" : This is the Miller-Madow asymptotic bias corrected empirical estimator.
#' \item "shrink" : This is a shrinkage estimate of the entropy of a Dirichlet probability distribution.
#' \item "sg" : This is the Schurmann-Grassberger estimate of the entropy of a Dirichlet probability distribution.
#' }
#' Type of Normalization:
#' \itemize{
#' \item marginal MI = 2*( Hx + Hy - Hxy ) / ( Hx + Hy ) 
#' \item joint     MI = 2*( Hx + Hy - Hxy ) / ( Hxy ) 
#' \item min.marginal	MI = ( Hx + Hy - Hxy ) / min(Hx,Hy) 
#' \item max.marginal MI = ( Hx + Hy - Hxy ) / max(Hx,Hy) 
#' \item min.conditional MI = ( Hx + Hy - Hxy ) / min(Hx.y,Hy.x) 
#' \item max.conditional MI = ( Hx + Hy - Hxy ) / max(Hx.y,Hy.x)
#' \item MI MI = ( Hx + Hy - Hxy )
#' }
#' @param X vector/factor denoting a random variable or a data.frame denoting a random vector where columns contain variables/features and rows contain outcomes/samples.
#' @param Y another random variable or random vector (vector/factor or data.frame).
#' @param method The name of the entropy estimator. The package implements four estimators : "emp", "mm", "shrink", "sg" (default:"emp") - see details. These estimators require discrete data values - see \code{\link{discretize}}.
#' @param type method of normalization. Default is "NULL" and the Mutual Information is calculated as MI = Hx+Hy-Hxy. Other methods include "MI", "marginal", "joint", "min.marginal", "max.marginal", "min.conditional", "max.conditional". See details below.
#' @return \code{nMI} returns the normalised mutual information I(X;Y) in nats.
#' 
#' @author Zhilong JIA
#' @seealso \code{\link{mulinformation}}, \code{\link[HDMD]{NMI}}
#' @references
#' {Meyer,  P. E.  (2008). Information-Theoretic Variable Selection and Network Inference from Microarray Data. PhD thesis of the Universite Libre de Bruxelles.}
#' 
#' {Cover, T. M. and Thomas, J. A. (1990). Elements of Information Theory. John Wiley, New York.}
#' 
#' {\code{\link[HDMD]{NMI}}}
#' @examples
#' data(USArrests)
#' dat<-discretize(USArrests)
#' 
#' I <- nMI(dat, method= "emp", type="max.marginal")
#' I2<- nMI(dat[,1],dat[,2])


nMI <-function(X, Y=NULL, method="emp", type="max.marginal")
{
    res <- NULL 
    if (is.null(Y))
        if(is.atomic(X))
            {stop("supply both 'X' and 'Y' or a matrix-like 'X'")
        } else {
        var.id <- NULL
        if( is.matrix(X) )
            X<-data.frame(X)
        if(is.data.frame(X)) 
            var.id <- names(X) 
        else stop("supply a matrix-like argument")
        
        X <- data.matrix(X)
        n <- NCOL(X)
        N <- NROW(X)
        Z<-na.omit(X)
        if( !(all(Z==round(Z))))
            stop("This method requires discrete values")                      
        #if(n>32000)
        #stop("too many variables")
        if( method == "emp")
            choi<-0
        else if( method == "mm" )
            choi<-1
        else if( method == "sg" )
            choi<-2
        else if(method == "shrink")
            choi<-3
        else stop("unknown method")
        
        ###normalisation type
        typ <- switch(type,
               marginal = 0,
               joint = 1,
               min.marginal = 2,
               max.marginal = 3,
               min.conditional = 4,
               max.conditional = 5,
               MI = 6,
               stop("unknown normalisation type!"))
        
        res <- .Call( "buildNMIM",X,N,n, choi, typ, PACKAGE="infotheo")
        dim(res) <- c(n,n)
        res <- as.matrix(res)
        rownames(res) <- var.id
        colnames(res) <- var.id
    } else {
        U<-data.frame(Y,X)
        Hxy<-entropy(U, method) 
        Hx<-entropy(X, method)
        Hy<-entropy(Y, method)
        
        Hy.x = Hxy - Hx
        Hx.y = Hxy - Hy
        res = switch(type, 
                         marginal = ifelse(Hx + Hy > 0, 2 * (Hx + Hy - Hxy)/(Hx + Hy), 0), 
                         joint = ifelse(Hxy > 0,  2 * (Hx + Hy - Hxy)/(Hxy), 0), 
                         min.marginal = ifelse(min(Hx, Hy) > 0, (Hx + Hy - Hxy)/min(Hx, Hy), 0), 
                         min.conditional = ifelse(min(Hx.y, Hy.x) > 0, (Hx + Hy - Hxy)/min(Hx.y, Hy.x), 0), 
                         max.conditional = ifelse(max(Hx.y, Hy.x) > 0, (Hx + Hy - Hxy)/max(Hx.y, Hy.x), 0), 
                         MI = Hx + Hy - Hxy,
                         #max.marginal
                         ifelse(max(Hx, Hy) > 0, (Hx + Hy - Hxy)/max(Hx, Hy), 0))
    }
    res
}