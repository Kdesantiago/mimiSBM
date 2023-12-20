#' Plot adjacency matrices
#' 
#' A function to plot each adjacency matrices defined by the thrid dimension of an array, and plot the sum of all theses matrices.
#' 
#' @param A an array with dim=c(N,N,V).
#'
#' @return None
#' @export
#' @importFrom graphics image
plot_adjacency<- function(A){
  #------------ Objectif ------------
  # Permet de visualiser les différentes matrices d'adjacences du tensor A

  #------------ Variable ------------
  # A : array(.,.,V)
  V <- dim(A)[3]

  #------------ Visualisation ------------
  for(v in 1:V){
    image(A[,,v],axes=F)
  }
  image(apply(A,c(1,2),sum))
}
