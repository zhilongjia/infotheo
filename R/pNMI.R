#' paralleled normalised mutual information computation
#'
#' Unsupervized Data Discretization and paralleled normalised mutual information computation. 
#'
#' method of the entropy estimator:
#' \itemize{
#' \item "emp" This estimator computes the entropy of the empirical probability distribution.
#' \item "mm" This is the Miller-Madow asymptotic bias corrected empirical estimator.
#' \item "shrink" This is a shrinkage estimate of the entropy of a Dirichlet probability distribution.
#' \item "sg" This is the Schurmann-Grassberger estimate of the entropy of a Dirichlet probability distribution.
#' }
#' type of Normalization
#' \itemize{
#' \item marginal MI = 2*( Hx + Hy - Hxy ) / ( Hx + Hy ) 
#' \item joint    MI = 2*( Hx + Hy - Hxy ) / ( Hxy ) 
#' \item min.marginal	MI = ( Hx + Hy - Hxy ) / min(Hx,Hy) 
#' \item max.marginal MI = ( Hx + Hy - Hxy ) / max(Hx,Hy) 
#' \item min.conditional MI = ( Hx + Hy - Hxy ) / min(Hx.y,Hy.x) 
#' \item max.conditional MI = ( Hx + Hy - Hxy ) / max(Hx.y,Hy.x)
#' \item MI = no Normalization
#' }
#'
#' @param X, vector/factor denoting a random variable or a data.frame denoting a 
#' random vector where \strong{column contain variables/features and row contain outcomes/samples.}
#' @param ncore number of cores
#' @param method The name of the entropy estimator. The package implements four 
#' estimators : "emp", "mm", "shrink", "sg" (default:"emp") - see details.
#' @param type method of normalization. Default is "max.marginal". Other methods 
#' include "marginal", "joint", "min.marginal", "max.marginal", "min.conditional", 
#' "MI", "max.conditional". See details below.
#' @param disc The name of the discretization method to be used :"equalfreq", 
#' "equalwidth" or "globalequalwidth" (default : "equalfreq") - see \code{discretize}.
#' @param nbins Integer specifying the number of bins to be used for the 
#' discretization. By default the number of bins is set to (N)^(1/3) where N is 
#' the number of samples.
#' 
#' 
#' @return a NMI matrix
#' @examples
#' data(USArrests)
#' pNMI(USArrests, 4)


pNMI <-function(X, ncore, method="emp", type="max.marginal", disc = "equalfreq", 
                nbins = NROW(X)^(1/3), verbose=FALSE){
    
    registerDoMC(ncore)
    
    #Unsupervized Data Discretization
    if (verbose) {print("discretize")}
    X <- foreach(i=1:ncol(X), .combine = cbind) %dopar% {
            discretize(X[,i], disc, nbins)
        }
    #X <- discretize(X, disc, nbins)
    if (verbose) {print("discretize done")}
    #registerDoMC(ncore)
    if (verbose) {print(paste("getDoParWorkers:", getDoParWorkers()))}
    NMI <- foreach(i=1:ncol(X), .combine=cbind) %:%
        foreach(j=1:ncol(X), .combine=c) %dopar% {
            Hxy <- entropy(data.frame(X[,c(i,j)]), method)
            Hx <- entropy(X[,i], method)
            Hy <- entropy(X[,j], method)
            Hy.x = Hxy - Hx
            Hx.y = Hxy - Hy
            NMI.val = switch(type, 
                             marginal = ifelse(Hx + Hy > 0, 2 * (Hx + Hy - Hxy)/(Hx + Hy), 0), 
                             joint = ifelse(Hxy > 0,  2 * (Hx + Hy - Hxy)/(Hxy), 0), 
                             min.marginal = ifelse(min(Hx, Hy) > 0, (Hx + Hy - Hxy)/min(Hx, Hy), 0), 
                             min.conditional = ifelse(min(Hx.y, Hy.x) > 0, (Hx + Hy - Hxy)/min(Hx.y, Hy.x), 0), 
                             max.conditional = ifelse(max(Hx.y, Hy.x) > 0, (Hx + Hy - Hxy)/max(Hx.y, Hy.x), 0), 
                             MI = Hx + Hy - Hxy,
                             #max.marginal
                             ifelse(max(Hx, Hy) > 0, (Hx + Hy - Hxy)/max(Hx, Hy), 0))
        }

    colnames(NMI) = rownames(NMI) = colnames(X)
    return (NMI)
}