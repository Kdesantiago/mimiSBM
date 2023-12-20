#---------------- Initialisation ----------------



#' SBM on each layer
#'
#' @param A an array of dim=c(N,N,V) 
#' @param silent Boolean for verbose
#' @param ncores the number of cores used to parallelize the calculations of the various SBMs
#'
#' @return a list containing the parameters of each SBM applied to each view
#' @importFrom blockmodels BM_bernoulli
#' @export
fit_SBM_per_layer <- function (A,silent=FALSE,ncores=2){
  M <- dim(A)[3]
  clustering <- vector("list", M)
  listTheta <- vector("list", M)
  for (m in 1:M){
    if (!silent){
      message("initialization of network", m)
    }
    myModel <- BM_bernoulli("SBM_sym", A[,,m], verbosity=0, plotting='',ncores=ncores)
    myModel$estimate()
    selectedModel <- which.max(myModel$ICL) + 1
    
    theta <- list(
      pi = colMeans(myModel$memberships[[selectedModel]]$Z),
      gamma = myModel$model_parameters[[selectedModel]]$pi,
      Z =  max.col(myModel$memberships[[selectedModel]]$Z),
      tau =  myModel$memberships[[selectedModel]]$Z)
    # permut <- degreeSort(theta, outTheta=FALSE, outPerm=TRUE)
    
    clustering[[m]] <- myModel$memberships[[selectedModel]]$Z
    #permutLabels(myModel$memberships[[selectedModel]]$map()$C,permut)
    
    # avoid empty blocks
    if(length(unique(clustering[[m]])) < max(unique(clustering[[m]]))){
      theta <- list(
        pi = colMeans(myModel$memberships[[selectedModel-1]]$Z),
        gamma = myModel$model_parameters[[selectedModel-1]]$pi,
        Z =  max.col(myModel$memberships[[selectedModel-1]]$Z),
        tau = myModel$memberships[[selectedModel-1]]$map()$C)
      #permut <- degreeSort(theta, outTheta=FALSE, outPerm=TRUE)
      
      clustering[[m]] <- myModel$memberships[[selectedModel-1]]$map()$C
      #permutLabels(myModel$memberships[[selectedModel-1]]$map()$C,  permut)
    }
    
    
    listTheta[[m]] <- theta
  }
  
  return(listTheta)
}

#' SBM on each layer - parallelized
#'
#' @param A an array of dim=c(N,N,V) 
#' @param nbCores the number of cores used to parallelize the calculations of the various SBMs
#'
#' @return a list containing the parameters of each SBM applied to each view
#' @importFrom parallel mclapply
#' @export
fit_SBM_per_layer_parallel <- function (A, nbCores=2){
  M <- dim(A)[3]
  nbCores <- min(c(ceiling(M/10), nbCores))
  Mpart <- ceiling(M/nbCores)
  init.values <- vector('list', nbCores)
  ind <- 1
  for (k in 1:nbCores){
    indEnd <- if (k<nbCores) ind + Mpart-1 else M
    init.values[[k]] <- A[,,ind:indEnd]
    ind <- ind + Mpart
  }
  
  res_parallel <- mclapply(init.values,
                           function(el){
                             fit_SBM_per_layer(el,silent=TRUE, ncores=nbCores)
                           },
                           mc.cores = nbCores
  )
  
  listTheta <- res_parallel[[1]]
  if(nbCores>1){
    for (k in 2:nbCores){
      listTheta <- c(listTheta, res_parallel[[k]])
    }
  }
  res <- listTheta
  
  return(res)
}


#' Initialization of mimiSBM parameters
#'
#' @param A an array of dim=c(N,N,V) 
#' @param K Number of clusters
#' @param Q Number of components
#' @param beta_0 hyperparameters for beta
#' @param theta_0 hyperparameters for theta
#' @param eta_0 hyperparameters for eta
#' @param xi_0 hyperparameters for xi
#' @param type_init select the type of initialization type_init=c("SBM","Kmeans","random")
#' @param nbCores the number of cores used to parallelize the calculations of the various SBMs 
#'
#' @return a list `params` updated
#' @importFrom stats kmeans runif
initialisation_params_bayesian <-function(A,K,Q,beta_0=rep(1/2,K),theta_0=rep(1/2,Q),eta_0=array(rep(1/2,K*K*Q),c(K,K,Q)),xi_0=array(rep(1/2,K*K*Q),c(K,K,Q)),type_init="SBM",nbCores=2){
  N = nrow(A)
  V = dim(A)[3]
  params <- list()
  
  params$beta_0 = beta_0
  params$theta_0 = theta_0
  params$eta_0 = eta_0
  params$xi_0 = xi_0
  
  if(type_init=="SBM"){
    if(nbCores>1) {
      simpleSBM <- fit_SBM_per_layer_parallel(A,nbCores = nbCores)
    } else {
      simpleSBM <- fit_SBM_per_layer(A,ncores = nbCores)
    }
    
    #----- Init vues
    mat <- c()
    for(v in 1:V){
      val <- sapply(seq(1:N),function(i){sapply(seq(1:N), function(j){ simpleSBM[[v]]$gamma[simpleSBM[[v]]$Z[i],simpleSBM[[v]]$Z[j]] } )})
      mat <- rbind(mat,as.vector(val))
    }
    params$u <- one_hot_errormachine(tryCatch({ kmeans(mat,centers=Q,nstart = 50)$cluster},
                                              error=function(cond){
                                                #message("Kmeans with small noise")
                                                mat2 <- mat + matrix(runif(n = nrow(mat) * ncol(mat),min=0,max=.1e-8),nrow(mat),ncol(mat))
                                                km= kmeans(mat2,centers=Q,nstart = 50)$cluster
                                                return(km)
                                              }))
    
    #----- Init clust individus
    
    mat <- c()
    for(v in (1:V)){
      mat <-cbind(mat,simpleSBM[[v]]$tau)
    }
    params$tau <- one_hot_errormachine(tryCatch({ kmeans(mat,centers=K,nstart = 50)$cluster},
                                                error=function(cond){
                                                  #message("Kmeans with small noise")
                                                  mat2 <- mat + matrix(runif(n = nrow(mat) * ncol(mat),min=0,max=.1e-8),nrow(mat),ncol(mat))
                                                  km= kmeans(mat2,centers=K,nstart = 50)$cluster
                                                  return(km)
                                                }))
  } else if(type_init=="Kmeans"){
    A_tmp <- apply(A,c(1,2),sum)
    params$tau <-tryCatch( {tmp = kmeans(A_tmp,centers = K,nstart = 50);
    one_hot_errormachine(tmp$cluster)},
    error=function(cond) {
      #message(paste("Bug KMeans, initialization with random parameters"))
      tmp = array(runif(N*K,0,1),dim=c(N,K)) ; tmp <- tmp/apply(tmp,1,sum)
      return(tmp)
    })
    
    params$u <- array(runif(N*K,0,1),dim=c(V,Q)) ; params$u <- params$u/apply(params$u,1,sum)
    
    
  } else {
    params$tau <- array(runif(N*K,0,1),dim=c(N,K)) ; params$tau <- params$tau/apply(params$tau,1,sum)
    params$u <- array(runif(N*K,0,1),dim=c(V,Q)) ; params$u <- params$u/apply(params$u,1,sum)
  }
  params <- update_beta_bayesian(params)
  params <- update_theta_bayesian(params)
  params <- update_eta_bayesian(A,params)
  params <- update_xi_bayesian(A,params)
  
  
  return(params)
}