
entropy <- function( X, estimator = "emp")
{
      X<-data.frame(X)
	  X<-data.matrix(X)
	  n <- NCOL(X)
      N <- NROW(X)
	  Z<-na.omit(X)
	  if( !(all(Z==round(Z))))
	      stop("This estimator requires discrete values")                      
	if(n>32000)
		stop("too many variables")
    res <- NULL 

    if( estimator == "emp")
		choi<-0
	else if( estimator == "mm" )
		choi<-1
	else if( estimator == "sg" )
		choi<-2
	else if(estimator == "shrink")
		choi<-3
	else stop("unknown estimator")
	res <- .Call( "entropyR",X,N,n, choi,
                        DUP=FALSE,PACKAGE="infotheo")
	res
}

multiinformation <- function( X, estimator = "emp")
{
	  X<-data.frame(X)
      X<-data.matrix(X)
	  n <- NCOL(X)
      N <- NROW(X)
	  Z<-na.omit(X)
	  if( !(all(Z==round(Z))))
	      stop("This estimator requires discrete values")                      
	if(n>32000)
		stop("too many variables")
    res <- NULL 

     if( estimator == "emp")
		choi<-0
	else if( estimator == "mm" )
		choi<-1
	else if( estimator == "sg" )
		choi<-2
	else if(estimator == "shrink")
		choi<-3
	else stop("unknown estimator")
	res <- .Call( "multiinformationR",X,N,n, choi,
                        DUP=FALSE,PACKAGE="infotheo")
	res
}

interinformation <- function( X, estimator = "emp")
{
	  X<-data.frame(X)
	  X<-data.matrix(X)
      n <- NCOL(X)
      N <- NROW(X)
	  
	 Z<-na.omit(X)
	  if( !(all(Z==round(Z))))
	      stop("This estimator requires discrete values")                      
	if(n>500)
		stop("too many variables")
    res <- NULL 

    if( estimator == "emp")
		choi<-0
	else if( estimator == "mm" )
		choi<-1
	else if( estimator == "sg" )
		choi<-2
	else if(estimator == "shrink")
		choi<-3
	else stop("unknown estimator")
	res <- .Call( "interactionR",X,N,n, choi,
                        DUP=FALSE,PACKAGE="infotheo")
	res
}

#compute H(X|Y)
condentropy<-function(X, Y=NULL, estimator="emp")
{
if(is.null(Y))
   Hres<-entropy(X, estimator)
else
   {
   Hyx<-entropy(data.frame(Y,X), estimator)
   Hy<-entropy(Y, estimator)
   Hres<-Hyx-Hy
   }
Hres
}

#compute I(X;Y)
mutinformation<-function(X, Y=NULL, estimator="emp")
{
	res <- NULL 
	if (is.null(Y))
		if(is.atomic(X))
			stop("supply both 'X' and 'Y' or a matrix-like 'X'")
		else {
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
				stop("This estimator requires discrete values")                      
			if(n>32000)
				stop("too many variables")
			if( estimator == "emp")
				choi<-0
			else if( estimator == "mm" )
				choi<-1
			else if( estimator == "sg" )
				choi<-2
			else if(estimator == "shrink")
				choi<-3
			else stop("unknown estimator")
			
			res <- .Call( "buildMIM",X,N,n, choi,
                        DUP=FALSE,PACKAGE="infotheo")
			dim(res) <- c(n,n)
			res <- as.data.frame(res)
			names(res) <- var.id
			row.names(res) <- var.id
			as.matrix(res)
		}
	else {
		U<-data.frame(Y,X)
		Hyx<-entropy(U, estimator) 
		Hx<-entropy(X, estimator)
		Hy<-entropy(Y, estimator)
		res<-Hx+Hy-Hyx
		if(res < 0)
			res<-0
	}
	res
}

#compute I(X;Y|S)
condinformation<-function(X,Y,S=NULL, estimator="emp")
{
if(is.null(S))
   Ires<-mutinformation(X,Y, estimator)
else
   {
   U<-data.frame(S,X,Y)
   Hysx<-entropy(U,estimator)
   Hsx<-entropy(U[,c(1,2)],estimator)
   Hys<-entropy(U[,c(1,3)],estimator)
   Hs<-entropy(S,estimator)
   Ires<- Hys - Hs - Hysx + Hsx
   }
Ires
}

natstobits<-function(H)
{ H*1.442695 }