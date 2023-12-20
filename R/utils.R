#' One Hot Encoding with Error machine
#'
#' @param Z a vector of size N, where Z\[i\] value indicate the cluster membership of observation i. 
#' @param size optional parameter, indicating the number of classes (avoid some empty class problems).
#'
#' @return Z a matrix N x K One-Hot-Encoded by rows, where K is the number of clusters.
#' @export
#' @examples
#' Z <- sample(1:4,10,replace=TRUE)
#' Z_OHE <- one_hot_errormachine(Z)
#' print(Z_OHE)
one_hot_errormachine <- function(Z,size=NULL){
  #------------ Objectif ------------
    # Prend en vecteur de classfication en entrée pour le One-Hot-Encode
    # en prenant en compte l'erreur machine possible.

  #------------ Variables ------------
  n <- length(Z)
  if(is.null(size)) K <- length(unique(Z)) else K = size

  #------------ One Hot Encoding + erreur machine ------------
  mat <- matrix(.Machine$double.xmin,n,K)
  mat[cbind(1:n,Z)] <- 1-.Machine$double.xmin

  #------------ Output ------------
  return(mat)
}

#' Diagonal coefficient to 0 on each slice given the 3rd dimension.
#'
#' @param A a array of dimension dim=c(N,N,V)
#'
#' @return A with 0 on each diagonal given the 3rd dimension.
#' @export
diag_nulle<- function(A){
  # Permet de récupérer un tenseur avec un coef 0 sur la diagonale des deux premières dimensions
  V <- dim(A)[3]
  for(v in 1:V) {diag(A[,,v]) <- 0 }
  return(A)
}


#'  Upper triangular Matrix/Array
#'
#' @param A a array or a squared matrix
#' @param transp boolean, indicate if we need a transposition or not.
#' @param diag boolean,  if True, diagonal is not used.
#'
#' @return a array or a squared matrix, with only upper-triangular coefficients with non-zero values
#' @export
trig_sup<-function(A,transp=F,diag=T){
  #------------ Objectif ------------
  # Permet de récupérer la matrice - où le tenseur - triangulaire supérieure

  #------------ Variables ------------
  # A : Matrice(.,.) ou array(.,.,V)

  #------------ Récupération triangulaire supérieure ------------
  if(length(dim(A))==3){
    V <- dim(A)[3]
    for(v in 1:V){
      tmp <- A[,,v]
      if(transp){tmp = t(tmp)}
      tmp[lower.tri(tmp,diag=diag)] <- 0
      A[,,v]<- tmp
    }
  } else {

  if(transp){A = t(A)}
  A[lower.tri(A,diag=T)] <- 0
}
  #------------ Output ------------
  return(A)
}

#' Transposition of an array
#'
#' @param A a array of dim= c(.,.,V)
#'
#' @return A_transposed, the transposed array according the third dimension
#' @export
transpo<- function(A){
  #------------ Objectif ------------
  # Permet de transposer un tenseur sur les deux premières dimensions

  #------------ transposition ------------
  V <- dim(A)[3]
  for(v in 1:V) {A[,,v] <- t(A[,,v]) }
  #------------ Output ------------
  return(A)
}

#' Sort the clustering matrix
#'
#' @param Z a matrix N x K, with probabilities to belong of a cluster in rows for each observation.
#'
#' @return a sorted matrix 
#' @export
sort_Z <- function(Z){
  #------------ Objectif ------------
  # Permet de réordonner un vecteur de labels (one-hot encoded ou non)

  #------------ Variables ------------
  # Z : vecteur de labels
  K = ncol(Z)

  #------------ Réordonnement ------------
  for(k in K:1){
    Z <- Z[order(Z[,k]),] # Réordonne les données
  }
  #------------ Output ------------
  return(Z)
}

#' Clustering Matrix : One hot encoding
#'
#' @param Z a matrix N x K, with probabilities to belong of a cluster in rows for each observation.
#'
#' @return  Z a matrix N x K One-Hot-Encoded by rows, where K is the number of clusters.
#' @export
#' @examples
#' Z <- matrix(rnorm(12),3,4)
#' Z_cem <- CEM(Z)
#' print(Z_cem)
CEM <- function(Z){
  #------------ Objectif ------------
  # transforme une matrice de label one-hot encoded en vecteur de labels

  #------------ Variable ------------
  n <- nrow(Z)

  #------------ Classfication EM ------------
  for(i in 1:n){
    tmp <- which.max(Z[i,])
    Z[i,tmp] = 1
    Z[i,-tmp] = 0
  }

  #------------ Output ------------
  return(Z)
}

#' log softmax of matrices (by row)
#'
#' @param log_X a matrix of log(X)
#' 
#' @return X with log_softmax function applied on each row
#' @importFrom stats rnorm
#' @export
#' @examples 
#' set.seed(42)
#' X <- matrix(rnorm(15,mean=5),5,3)
#' log_X <- log(X)
#' X_softmax <-  log_Softmax(X)
log_Softmax <- function(log_X){
  if(!is.matrix(log_X)){log_X <- as.matrix(log_X)}
  K <- ncol(log_X)

  log_X <- log_X - apply(log_X,1,max)

  ## Now going back to exponential with the same normalization
  X <- exp(log_X) #(matrix(1,n,1) %*% pi) * exp(logX)
  X <- pmin(X,.Machine$double.xmax)
  X <- pmax(X,.Machine$double.xmin)
  X <- X / (rowSums(X) %*% matrix(1,1,K))
  X <- pmin(X,1-.Machine$double.xmin)
  X <- pmax(X,.Machine$double.xmin)

  return(X)
}

#' Calculation of Log multinomial Beta value.
#'
#' 
#' @param x a vector
#'
#' @return sum(lgamma(x\[j\])) - lgamma(sum(x))
#' @export 
multinomial_lbeta_function <- function(x){sum(lgamma(x))-lgamma(sum(x))}
