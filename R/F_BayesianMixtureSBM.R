require(blockmodels)
require(parallel)

#---------------- Update parameters ----------------

update_tau_bayesian<-function(A,params,eps_conv=1e-4){

  tau= params$tau
  u = params$u
  beta_k = params$beta_k
  eta = params$eta
  xi = params$xi

  N = dim(A)[1]
  K = ncol(tau)
  Q = ncol(u)
  V = nrow(u)
  log_tau <- matrix(0,N,K)

  A <- diag_nulle(A)
  old_tau <- tau + 1

  etp <- 0

  while( (sum(abs(old_tau - tau)) > eps_conv) && etp<50 ){
    etp <- etp+1
    old_tau <- tau

    for(i in 1:N){
      tau_tmp = tau
      tau_tmp[i,] = 0 # pour enlever le cas où i = j
      for(k in 1:K){
        for(q in 1:Q){
          log_tau[i,k] = log_tau[i,k] +  t(u[,q]) %*%(t(A[i,,]) %*% tau) %*%  ( digamma(eta[k,,q]) - digamma(xi[k,,q])) #Partie Aijv de la somme
        }
        log_tau[i,k] = log_tau[i,k] + sum( tau_tmp %*% (digamma(xi[k,,]) - digamma(eta[k,,] + xi[k,,])) %*% t(u)) # Deuxième partie (sans Aijv) de la somme

        log_tau[i,k] = log_tau[i,k] + digamma(beta_k[k]) - digamma(sum(beta_k)) # Partie esp(pi_k)
      }
    }

    tau <- log_Softmax(log_tau)

  }
  params$tau = tau
  return(params)
}

update_u_bayesian<-function(A,params){

  tau= params$tau
  theta = params$theta
  eta = params$eta
  xi = params$xi

  V = dim(A)[3]
  Q = length(theta)
  K = ncol(tau)
  log_u <- matrix(0,V,Q)

  A <- diag_nulle(A)
  A_trig_sup <- trig_sup(A)


  for(q in 1:Q){
    for(v in 1:V){
      for(k in 1:K){
        for(l in k:K){
          if(l == k){
            log_u[v,q] = log_u[v,q] + tau[,k] %*%  ( A_trig_sup[,,v]*( digamma(eta[k,l,q]) - digamma(xi[k,l,q])) + digamma(xi[k,l,q]) - digamma(eta[k,l,q] + xi[k,l,q]) ) %*% tau[,l]
          } else {
            log_u[v,q] = log_u[v,q] + tau[,k] %*%  (A[,,v]*(digamma(eta[k,l,q]) - digamma(xi[k,l,q])) + digamma(xi[k,l,q]) - digamma(eta[k,l,q] + xi[k,l,q]) ) %*% tau[,l]
          }
        }
      }
      log_u[v,q] = log_u[v,q] + digamma(theta[q]) - digamma(sum(theta))
    }
  }

  u <- log_Softmax(log_u)

  params$u <- u
  return(params)
}

update_beta_bayesian<-function(params){

  tau= params$tau
  beta_0 = params$beta_0

  beta_k <- beta_0 + apply(tau,2,sum)
  params$beta_k <- beta_k
  return(params)
}

update_theta_bayesian<-function(params){

  u= params$u
  theta_0 = params$theta_0

  theta <- theta_0 + apply(u,2,sum)
  params$theta <- theta
  return(params)
}

update_eta_bayesian<-function(A,params){

  eta_0 = params$eta_0
  tau= params$tau
  u= params$u

  K = ncol(tau)
  Q = ncol(u)
  eta <- array(0,c(K,K,Q))
  A = diag_nulle(A)
  A_trig_sup = trig_sup(A)

  for(q in 1:Q){
    for(k in 1:K){
      for(l in k:K){
        if(l == k){
          eta[k,l,q] = tau[,k] %*%  apply(A_trig_sup,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
        } else {
          eta[k,l,q] = tau[,k] %*%  apply(A,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
          eta[l,k,q] = eta[k,l,q]
        }
      }
    }
    eta[k,l,q] = eta[k,l,q] + eta_0[k,l,q]
  }
  params$eta <- pmax(eta,.Machine$double.eps)#eta
  return(params)
}

update_xi_bayesian<-function(A,params){

  xi_0 = params$xi_0
  tau= params$tau
  u= params$u

  K = ncol(tau)
  Q = ncol(u)
  xi <- array(0,c(K,K,Q))

  A = 1-A # Car xi dépend de 1 - Aijv
  A = diag_nulle(A)
  A_trig_sup = trig_sup(A)

  for(q in 1:Q){
    for(k in 1:K){
      for(l in k:K){
        if(l == k){
          xi[k,l,q] = t(tau[,k]) %*%  apply(A_trig_sup,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
        } else {
          xi[k,l,q] = t(tau[,k]) %*%  apply(A,c(1,2),function(x){x %*% u[,q]}) %*% tau[,l]
          xi[l,k,q] = xi[k,l,q]
        }
      }
    }
    xi[k,l,q] = xi[k,l,q] + xi[k,l,q]
  }
  params$xi <-  pmax(xi,.Machine$double.eps) #xi
  return(params)
}

multinomial_lbeta_function <- function(x){sum(lgamma(x))-lgamma(sum(x))}

#---------------- Loss ----------------
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


#---------------- Variational Bayes EM ----------------

VBEM_step<-function(A,params,alternate=T,eps_conv=1e-3){

  if(alternate){
    params <- update_u_bayesian(A,params)
    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)

    params <- update_tau_bayesian(A,params,eps_conv)
    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)


  } else{
    params <- update_u_bayesian(A,params)
    params <- update_tau_bayesian(A,params,eps_conv)
    params <- update_beta_bayesian(params)
    params <- update_theta_bayesian(params)
    params <- update_eta_bayesian(A,params)
    params <- update_xi_bayesian(A,params)
  }


  return(params)
}

#---------------- Initialisation ----------------




#' import
fit_SBM_per_layer <- function (A,silent=FALSE,ncores=detectCores()){
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

fit_SBM_per_layer_parallel <- function (A, nbCores=detectCores()){
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


initialisation_params_bayesian <-function(A,K,Q,beta_0=rep(1/2,K),theta_0=rep(1/2,Q),eta_0=array(rep(1/2,K*K*Q),c(K,K,Q)),xi_0=array(rep(1/2,K*K*Q),c(K,K,Q)),type_init="SBM"){
  N = nrow(A)
  V = dim(A)[3]
  params <- list()

  params$beta_0 = beta_0
  params$theta_0 = theta_0
  params$eta_0 = eta_0
  params$xi_0 = xi_0

  if(type_init=="SBM"){
    if(detectCores()>1) {
      simpleSBM <- fit_SBM_per_layer_parallel(A)
    } else {
      simpleSBM <- fit_SBM_per_layer(A)
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



#---------------- Model ----------------

BayesianMixture_SBM_model <-function(A,K,Q,beta_0=rep(1/2,K),theta_0=rep(1/2,Q),eta_0=array(rep(1/2,K*K*Q),c(K,K,Q)),xi_0=array(rep(1/2,K*K*Q),c(K,K,Q)),tol=1e-3,iter_max=10,n_init = 1,alternate=T, Verbose=T,eps_conv=1e-4,type_init="SBM"){


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
    params <- initialisation_params_bayesian(A,K=K,Q=Q,beta_0,theta_0,eta_0,xi_0,type_init)
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
      if(Verbose){print(paste0("__ Intération n°", n_iter," __"))}

      params <- VBEM_step(A,params,alternate,eps_conv)
      elbo <- Loss_BayesianMSBM(params)


      if(Verbose){print(paste0("L'Evidence Lower BOund vaut :",round(elbo,2)))}
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





BayesianMixture_SBM <-function(A,Kset,Qset,criterion = "ILVB",tol=1e-3,iter_max=10,n_init = 1,alternate=T, Verbose=F,eps_conv=1e-4,type_init="SBM", nbCores = detectCores()-1){
  val_crit_best = -Inf
  val_crit = -Inf
  model_best = NULL

  N = dim(A)[1]
  V = dim(A)[3]
  for(q in Qset){
    for(k in Kset){
      #if(Verbose){print(paste0("___K : ",k, " et Q : ",q," ___"))}
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
                                 print(el)
                                 BayesianMixture_SBM_model(A=A,K=k,Q=q,iter_max=iter_max,tol=tol,n_init=length(el),alternate=alternate, Verbose=Verbose,eps_conv=eps_conv,type_init=type_init)
                               },
                               mc.cores = nbCores
      )

      idx = which.max(sapply(1:nbCores, function(i){res_parallel[[i]]$elbo})) # A modifier
      model = res_parallel[[idx]]
      model$Q <- q
      model$K <- k
      #model = BayesianMixture_SBM_model(A=A,K=k,Q=q,iter_max=iter_max,tol=tol,n_init=n_init,alternate=alternate, Verbose=Verbose,eps_conv=eps_conv,type_init=type_init)$best

      message("Components verification")
      warn = (FALSE %in% sapply(1:k, function(i){ i %in% max.col(model$tau)})) || (FALSE %in% sapply(1:q, function(i){ i %in% max.col(model$u)}) )
      if(warn){
        warning("probleme dans les hyperparamètres K ou Q choisi")
        #next

        if(FALSE %in% sapply(1:q, function(i){ i %in% max.col(model$u)})){
          idx = which(!sapply(1:q, function(i){ i %in% max.col(model$u)}))
          model$u <- model$u[,-idx]
          model$u <- t(apply(model$u,1,function(i){i/sum(i)}))
          model$theta_0 <- model$theta_0[-idx]
          model$theta <-model$theta[-idx]
          model$eta_0 <- model$eta_0[,,-idx]
          model$eta <- model$eta[,,-idx]
          model$xi_0 <-  model$xi_0[,,-idx]
          model$xi <-  model$xi[,,-idx]
          model$Q <- model$Q - length(idx)
        }

        if(FALSE %in% sapply(1:k, function(i){ i %in% max.col(model$tau)})){
          idx = which(!sapply(1:k, function(i){ i %in% max.col(model$tau)}))
          model$tau <- model$tau[,-idx]
          model$tau <- t(apply(model$tau,1,function(i){i/sum(i)}))
          model$beta_0 <- model$beta_0[-idx]
          model$beta_k <-model$beta_k[-idx]
          model$eta_0 <- model$eta_0[-idx,-idx,]
          model$eta <- model$eta[-idx,-idx,]
          model$xi_0 <-  model$xi_0[-idx,-idx,]
          model$xi <-  model$xi[-idx,-idx,]
          model$K <- model$K - length(idx)
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
