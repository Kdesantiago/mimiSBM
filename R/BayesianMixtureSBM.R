#---------------- Loss ----------------
#' mimiSBM Evidence Lower BOund
#'
#' @param params  alist of parameters of the model
#'
#' @return Computation of the mimiSBM ELBO
#' @export
Loss_BayesianMSBM <- function(params){
  tau = params$tau
  u = params$u
  beta_k = pmax(params$beta_k,.Machine$double.xmin)
  theta = pmax(params$theta,.Machine$double.xmin)
  eta = params$eta
  xi = params$xi
  beta_0 = params$beta_0
  theta_0 = params$theta_0
  eta_0 = params$eta_0
  xi_0 <- params$xi_0

  K = ncol(tau)
  Q = ncol(u)

  res <- - sum(tau*log(tau)) - sum(u*log(u))

  res <- res + multinomial_lbeta_function(beta_k) - multinomial_lbeta_function(beta_0)


  res <- res + multinomial_lbeta_function(theta) - multinomial_lbeta_function(theta_0)

  tmp1 <- array(NA,c(K,K,Q))
  tmp2 <- array(NA,c(K,K,Q))
  for(k in 1:K){
    for(l in 1:K){
      for(q in 1:Q){
        tmp2[k,l,q] = beta(eta_0[k,l,q],xi_0[k,l,q])
        tmp1[k,l,q] = beta(eta[k,l,q],xi[k,l,q])
      }
    }
  }
  tmp2 <- pmax(tmp2,.Machine$double.xmin)
  tmp1 <- pmax(tmp1,.Machine$double.xmin)

  res <- res + sum(trig_sup(log(tmp1 / tmp2),diag = F))

  return(res)

}


#---------------- Model ----------------

#' mimiSBM model for fixed K and Q 
#'
#' @param A an array of dim=c(N,N,V) 
#' @param K Number of clusters
#' @param Q Number of components
#' @param beta_0 hyperparameters for beta
#' @param theta_0 hyperparameters for theta
#' @param eta_0 hyperparameters for eta
#' @param xi_0 hyperparameters for xi
#' @param tol convergence parameter on ELBO
#' @param iter_max maximal number of iteration of mimiSBM
#' @param n_init number of initialization of the mimi algorithm.
#' @param alternate boolean indicated if we put an M-step after each part of the E-step, after u optimization and after tau optimization. If not, we optimize u and tau and after the M-step is made. 
#' @param Verbose boolean for information on model fitting
#' @param eps_conv parameter of convergence for tau.
#' @param type_init select the type of initialization type_init=c("SBM","Kmeans","random")
#' @param nbCores the number of cores used to parallelize the calculations
#'
#' See the vignette for more details
#'
#' @return mimiSBM model with estimation of coefficients.
#' @export
#' 
BayesianMixture_SBM_model <-function(A,K,Q,beta_0=rep(1/2,K),theta_0=rep(1/2,Q),eta_0=array(rep(1/2,K*K*Q),c(K,K,Q)),xi_0=array(rep(1/2,K*K*Q),c(K,K,Q)),tol=1e-3,iter_max=10,n_init = 1,alternate=T, Verbose=T,eps_conv=1e-4,type_init="SBM",nbCores=2){


  #------------ Variables ------------
  N = nrow(A)
  V = dim(A)[3]
  output <- list()
  output$parametres <-list()
  output$n_iter <- rep(NA,n_init)
  output$elbo <- rep(NA,n_init)


  #------------ Boucle d'initialisation ------------
  for(init in 1:n_init)  {
    if(Verbose){print(paste0("------------ Run ",init," ------------"))}
    n_iter = 0
    params <- initialisation_params_bayesian(A,K=K,Q=Q,beta_0,theta_0,eta_0,xi_0,type_init,nbCores=nbCores)
    elbo = Loss_BayesianMSBM(params)


    elbo_old <- -Inf
    params.old <- params

    params_best <- params
    elbo_best <- -Inf

    #------------ Boucle VBEM - 1 run ------------
    while( (n_iter < iter_max) && (abs(elbo_old-elbo)>tol) ){
      elbo_old <- elbo
      params.old <- params

      n_iter = n_iter + 1
      if(Verbose){print(paste0("__ Interation n", n_iter," __"))}

      params <- VBEM_step(A,params,alternate,eps_conv)
      elbo <- Loss_BayesianMSBM(params)


      if(Verbose){print(paste0("Evidence Lower BOund is:",round(elbo,2)))}
      #if(elbo < elbo_old){ warning("l'ELBO diminue") ;}# params <- params.old}
      if(elbo>elbo_best){elbo_best <- elbo;params_best <- params}
    }

    output$parametres[[init]] <- params_best
    output$n_iter[init] <- n_iter
    output$elbo[init] <- elbo_best
  }


  #------------ Output ------------

  #output$best
  output$best <- output$parametres[[which.max(output$elbo)]]
  output$best$elbo <- output$elbo[which.max(output$elbo)]
  return(output$best)

}





#' Mixture of Multilayer Integrator SBM (mimiSBM)
#' 
#' mimiSBM is a model that allows both clustering of individuals and grouping of views by component.
#' This Bayesian model estimates the probability of individuals belonging to each cluster (cluster crossing all views) and the membership component for all views.
#' In addition, the connectivity tensor between classes, conditional on the components, is also estimated. 
#' @param A an array of dim=c(N,N,V) 
#' @param Kset Set of number of clusters
#' @param Qset Set of number of components
#' @param beta_0 hyperparameters for beta
#' @param theta_0 hyperparameters for theta
#' @param eta_0 hyperparameters for eta
#' @param xi_0 hyperparameters for xi
#' @param criterion  model selection criterion, criterion=c("ILVB","ICL_approx","ICL_variationnel","ICL_exact")
#' @param tol convergence parameter on ELBO
#' @param iter_max maximal number of iteration of mimiSBM
#' @param n_init number of initialization of the mimi algorithm.
#' @param alternate boolean indicated if we put an M-step after each part of the E-step, after u optimization and after tau optimization. If not, we optimize u and tau and after the M-step is made. 
#' @param Verbose boolean for information on model fitting
#' @param eps_conv parameter of convergence for tau.
#' @param type_init select the type of initialization type_init=c("SBM","Kmeans","random")
#'
#' @return The best model, conditionnally to the criterion, and its parameters.
#' @export
#'
#' @examples

mimiSBM <-function(A,Kset,Qset,beta_0=1/2,theta_0=1/2,eta_0=1/2,xi_0=1/2,criterion = "ILVB",tol=1e-3,iter_max=10,n_init = 1,alternate=T, Verbose=F,eps_conv=1e-4,type_init="SBM"){
  nbCores = 1
  val_crit_best = -Inf
  val_crit = -Inf
  model_best = NULL

  N = dim(A)[1]
  V = dim(A)[3]
  for(q in Qset){
    for(k in Kset){
      if(Verbose){print(paste0("___K : ",k, " et Q : ",q," ___"))}
      #print(paste0("___Q : ",q, " et K : ",k," ___"))

      #------- PARALLELISATION --------

      nbCores <- min(c(n_init, nbCores))
      Mpart <- ceiling(n_init/nbCores)
      init.values <- vector('list', nbCores)
      ind <- 1
      for (j in 1:nbCores){
        indEnd <- if (j<nbCores) ind + Mpart-1 else n_init
        init.values[[j]] <- ind:indEnd
        ind <- ind + Mpart
      }

      res_parallel <- mclapply(init.values,
                               function(el){
                                 #print(el)
                                 BayesianMixture_SBM_model(A=A,K=k,Q=q,beta_0=rep(beta_0,k),theta_0=rep(theta_0,q),eta_0=array(rep(eta_0,k*k*q),c(k,k,q)),xi_0=array(rep(xi_0,k*k*q),c(k,k,q)),iter_max=iter_max,tol=tol,n_init=length(el),alternate=alternate, Verbose=Verbose,eps_conv=eps_conv,type_init=type_init,nbCores = nbCores)
                               },
                               mc.cores = nbCores
      )

      idx = which.max(sapply(1:nbCores, function(i){res_parallel[[i]]$elbo})) # A modifier
      model = res_parallel[[idx]]
      model$Q <- q
      model$K <- k
      #model = BayesianMixture_SBM_model(A=A,K=k,Q=q,iter_max=iter_max,tol=tol,n_init=n_init,alternate=alternate, Verbose=Verbose,eps_conv=eps_conv,type_init=type_init)$best

      message("Components verification")
      Z_est = max.col(model$tau)
      W_est = max.col(model$u)
      warn = (FALSE %in% sapply(1:k, function(i){ i %in% Z_est})) || (FALSE %in% sapply(1:q, function(i){ i %in% W_est}) )
      if(warn){
        warning("probleme dans les hyperparametres K ou Q choisi")
        #next

        if(FALSE %in% sapply(1:q, function(i){ i %in% W_est})){
          idx = which(!sapply(1:q, function(i){ i %in% W_est}))
          model$Q <- model$Q - length(idx)
          if(model$Q == 1) {
            #print("ATTENTION Q = 1")
            model$u <- matrix(model$u[,-idx],V,1)
            model$theta_0 <- model$theta_0[-idx]
            model$theta <- model$theta[-idx]
            model$eta_0 <- array(model$eta_0[,,-idx],dim=c(model$K,model$K,1))
            model$eta <- array(model$eta[,,-idx],dim=c(model$K,model$K,1))
            model$xi_0 <-  array(model$xi_0[,,-idx],dim=c(model$K,model$K,1))
            model$xi <-  array(model$xi[,,-idx],dim=c(model$K,model$K,1))
          } else {
            model$u <- model$u[,-idx]
            model$u <- t(apply(model$u,1,function(i){i/sum(i)}))
            model$theta_0 <- model$theta_0[-idx]
            model$theta <-model$theta[-idx]
            model$eta_0 <- array(model$eta_0[-idx,-idx,],dim=c(model$K,model$K,model$Q))
            model$eta <- array(model$eta[-idx,-idx,],dim=c(model$K,model$K,model$Q))
            model$xi_0 <-  array(model$xi_0[-idx,-idx,],dim=c(model$K,model$K,model$Q))
            model$xi <-  array(model$xi[-idx,-idx,],dim=c(model$K,model$K,model$Q))
          }


        }

        if(FALSE %in% sapply(1:k, function(i){ i %in% Z_est})){
          idx = which(!sapply(1:k, function(i){ i %in% Z_est}))
          model$K <- model$K - length(idx)
          if(model$K == 1) {
            #print("ATTENTION K = 1")
            model$tau <- matrix(model$tau[,-idx],N,1)
            model$beta_0 <- model$beta_0[-idx]
            model$beta_k <-model$beta_k[-idx]
            model$eta_0 <- array(model$eta_0[-idx,-idx,],dim=c(1,1,model$Q))
            model$eta <- array(model$eta[-idx,-idx,],dim=c(1,1,model$Q))
            model$xi_0 <-  array(model$xi_0[-idx,-idx,],dim=c(1,1,model$Q))
            model$xi <-  array(model$xi[-idx,-idx,],dim=c(1,1,model$Q))
          } else {
            model$tau <- model$tau[,-idx]
            model$tau <- t(apply(model$tau,1,function(i){i/sum(i)}))
            model$beta_0 <- model$beta_0[-idx]
            model$beta_k <-model$beta_k[-idx]
            model$eta_0 <- array(model$eta_0[-idx,-idx,],dim=c(model$K,model$K,model$Q))
            model$eta <- array(model$eta[-idx,-idx,],dim=c(model$K,model$K,model$Q))
            model$xi_0 <-  array(model$xi_0[-idx,-idx,],dim=c(model$K,model$K,model$Q))
            model$xi <-  array(model$xi[-idx,-idx,],dim=c(model$K,model$K,model$Q))

          }


        }

      }


      if(criterion == "ICL_approx" ){
        val_crit = model$elbo - 1/2 * (k*(k+1)/2*q)*log(N*(N-1)/2*V) - 1/2 * (k-1) * log(N) - 1/2 * (q-1) * log(V)
      } else if(criterion == "ICL_variationnel"){
        val_crit = model$elbo + sum(model$tau*log(model$tau)) + sum(model$u*log(model$u))
      } else if(criterion == "ICL_exact"){

        params_tamp = list(tau = one_hot_errormachine(max.col(model$tau)),u = one_hot_errormachine(max.col(model$u)) )
        params_tamp$beta_0 = model$beta_0
        params_tamp$theta_0 = model$theta_0
        params_tamp$eta_0 = model$eta_0
        params_tamp$xi_0 <- model$xi_0
        params_tamp <- update_beta_bayesian(params_tamp)
        params_tamp <- update_theta_bayesian(params_tamp)
        params_tamp <- update_eta_bayesian(A,params_tamp)
        params_tamp <- update_xi_bayesian(A,params_tamp)
        val_crit = model$elbo + sum(model$tau*log(model$tau)) + sum(model$u*log(model$u))
      } else{ #"ILVB"
        val_crit = model$elbo
      }

      if(val_crit_best< val_crit){
        model_best = model
        val_crit_best = val_crit
      }

    }
  }

 model_best$criterion <- criterion
 model_best$val_criterion <- val_crit_best

 return(model_best)
}
