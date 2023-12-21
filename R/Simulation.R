#' Label Switching
#'
#'This function can be used to perturb a clustering vector in order to randomly associate certain individuals with another cluster.
#'
#' @param Z a clustering vector
#' @param p_out a probability of perturbation for the clustering
#'
#' @return a perturbed clustering vector
#' @importFrom stats runif
#' @export
#'
#' @examples
#' Z <- sample(1:4,100,replace=TRUE)
#' p = 0.1
#' Z_pert <- lab_switching(Z,p)
#' table("Initial clustering" = Z,"Perturbed clustering" = Z_pert)
lab_switching <- function(Z,p_out=0.1){
  N = length(Z)
  Z_tilde = Z
  val = runif(N)
  idx = which(val < p_out)
  Z_tilde[idx] <- sapply(idx, function(i) {
    set_clust <- unique(Z)
    set_clust <- setdiff(set_clust,Z[i])
    Z_tilde[i] = sample(set_clust,1)
  })

  return(unlist(Z_tilde))
}





#' Simulate data from the mimiSBM generative model.
#'
#' @param N Number of individuals.
#' @param V Number of views.
#' @param alpha_klq array of component-connection probability (K,K,Q).
#' @param pi_k Vector of proportions of individuals across clusters.
#' @param rho Vector of proportion of views across components.
#' @param sorted Boolean for simulation reordering (clusters and components membership).
#' @param p_switch probability of label-switching, if NULL no perturbation between true clustering and the connectivity of individuals.
#' @return list with the parameters of the simulation ($params), and
#'  the simulations ($simulation).
#' @importFrom stats rmultinom rbinom
#' @export
#---------------- Simulations basées sur le modèle ----------------
rMSBM <- function(N, V,alpha_klq,pi_k,rho,sorted=TRUE,p_switch=NULL){
  #------------ Objectif ------------
  # Simuler des données à partir du modèle génératif mimiSBM

  #------------ Variables ------------
  # N : Nomber d'individus
  # K : Nombre de  clusters
  # V : Nombre de vues
  # Q : Nombre de composantes du mélange de vue

  # alpha_klq : array(K,K,Q) Tenseur de probabilité de lien
  # pi_k : vecteur(K) Probabilité d'appartenir aux clusters
  # rho : vecteur(Q) Probabilité d'appartenir au mélange de vue Q
  #list(params = params, simulation=simulation)
  # Z : Variable latente individus -> appartenance aux clusters
  # W : Variable latente vues -> d'appartenance aux composantes de vues

  K = length(pi_k)
  Q = length(rho)
  A = rep(NA,N*N*V)
  attr(A,"dim") <- c(N,N,V)
  #------------ états cachés clustering ------------
  Z <- t(rmultinom(n=N, size=1, prob=pi_k))
  Z = apply(Z,1,which.max)
  if(sorted) {Z = Z[order(Z)]}

  #------------ mélange vues ------------
  W <- t(rmultinom(n=V, size=1, prob=rho))
  W <- apply(W,1,which.max)
  if(sorted) {W = W[order(W)]}

  #------------ Tenseur d'adjacence ------------

  for (v in 1:V){

    if(is.null(p_switch)) {
      Z_switch = Z
    } else {
      Z_switch = lab_switching(Z,p_switch)
    }

    for(i in 1:(N-1)){
      for(j in (i+1):N){
        A[i,j,v] <- rbinom(1,1,prob=alpha_klq[Z_switch[i],Z_switch[j],W[v]])
        A[j,i,v] <- A[i,j,v]
      }
    }
    diag(A[,,v])<- 0
  }

  #------------ output ------------

  params <- list(N=N,K=K,V=V,Q=Q,pi_k=pi_k,rho=rho,alpha_klq=alpha_klq)
  simulation <- list(A=A,Z=Z,W=W)
  output <- list(params = params, simulation=simulation)
  return(output)
}


#---------------- Simulations Recouvrement ----------------


#' Create a link between final clustering and clustering per view component.
#'
#' @param K_barre Number of clusters in the final clustering 
#' @param K Vector of size Q, indicate the number of clusters in each component. K\[q\] <= K_barre for all q
#' @return cluster : a list of link between final clustering and clustering per view component.
#' @export
partition_K_barre <-function(K_barre,K){
  #---------------------- Variables ----------------------------
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # K_barre : Nombre de vrais clusters
  ## Hypothèse : K[q] <= K_barre pour tout q

  if(any(K > K_barre)) {errorCondition("The number of clusters in at least 1 component is greater than K_barre")}
  Q = length(K)
  clusters <- list()

  #---------------------- Partition clusters ----------------------------
  ## On commence par associer 1 classe à chaque cluster de K[q]
  for(q in 1:Q){
    tmp <- list()
    idx_K_barre = 1:K_barre
    idx <- sample(idx_K_barre,replace = F,size=K[q])
    idx_K_barre = idx_K_barre[-which(idx_K_barre %in% idx)]
    for(l in 1:K[q]){tmp[[l]] = idx[l]}

    ## On fait un tirage aléatoire pour chaque classe qui n'a pas été associée.
    for(idx in idx_K_barre){
      idx_K <- sample(1:K[q],size=1)
      tmp[[idx_K]] <- c(tmp[[idx_K]],idx)
    }

    clusters[[q]] <- tmp
  }

  # ----------------------  Output ----------------------------

  return(clusters) # Cluster[[Q]] [[K[Q]]]

}


#' Create probality-component list for clustering per view component.
#' @param clusters list of link between final clustering and clustering per view component.
#' @param K_barre Number of clusters in the final clustering
#' @param K Vector of size Q, indicate the number of clusters in each component.
#' @return alpha :  probality-component list for clustering per view component.
#' @export
Mat_lien_alpha <- function(clusters,K_barre,K){

  #---------------------- Variables ----------------------------
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # K_barre : Nombre de vrais clusters
  # clusters : list indiquant les liens entre K_Barre et K ( Cluster[[Q]][[K[Q]]]  )

  Q = length(K)

  #---------------------- Simulation alpha ----------------------------
  alpha <- array(0.01,dim=c(K_barre,K_barre,Q))
  for(q in 1:Q){
    for(s in 1:K[q]){
      for(k in clusters[[q]][[s]]){
        for(l in clusters[[q]][[s]]){
          alpha[k,l,q] = 0.99
        }
      }
    }
  }
  return(alpha)
}



#' Simulation of mixture multilayer Stochastick block model
#'
#' This simulation process assumes that we have partial information on the clustering within each view component, and that the final clustering of individuals depends on a combination of the clustering on each of the views.
#' In addition, we take into account possible label-switching: we consider that an individual belongs with a certain probability to the wrong class, thus disturbing the adjacency matrices and making the simulation more real and complex.
#' 
#' See the vignette for more information.
#' @param N Number of observations
#' @param V Number of views
#' @param K Vector of size Q, indicate the number of clusters in each component.
#' @param pi_k Vector of proportions of observations across clusters.
#' @param rho Vector of proportion of views across components.
#' @param sorted Boolean for simulation reordering (clusters and components membership).
#' @param p_switch probability of label-switching, if NULL no perturbation between true clustering and the connectivity of individuals.
#'
#' @return list with the parameters of the simulation ($params), and the simulations ($simulation).
#' @export
rSMB_partition <- function(N, V, K, pi_k,rho,sorted=T,p_switch=NULL){

  #---------------------- Variables ----------------------------
  # N : Nombre d'observations à générer
  # V : Nombre de vues à générer
  # K : Nombre de clusters au sein de chaque composantes  (vecteur de taille Q)
  # K_barre : Nombre de vrais clusters
  # pi_k : Vecteur de probabilité d'appartenance aux clusters de K_barre
  # rho : Vecteur de probabilité d'appartenance aux clusters de W (mélange vue)
  # sorted : Boolean pour le réagencement des simulations.

  # clusters : list indiquant les liens entre K_Barre et K ( Cluster[[Q]][[K[Q]]]  )
  # alpha : array, probabilité de liaison entre 2 communautés. (dim = K_Barre,K_barre,Q)

  K_barre <- length(pi_k)

  #----------------------  K_Barre -> K ----------------------------

  clusters <- partition_K_barre(K_barre,K)

  # ----------------------  Alpha ----------------------------

  alpha <- Mat_lien_alpha(clusters,K_barre,K)

  # ----------------------  Z / W / A ----------------------------

  output <- rMSBM(N,V,alpha_klq = alpha,pi_k,rho,sorted=sorted, p_switch=p_switch)

  # ----------------------  Output ----------------------------

  output$params$clusters <- clusters
  output$params$K_barre <- K_barre
  output$params$K <- K
  return(output)
}




