#' Update of bayesian parameter tau
#'
#' @param A an array of dim=c(N,N,V) 
#' @param params list of parameters of the model
#' @param eps_conv parameter of convergence.
#'
#' @return params with tau updated
#' @export
#'
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

#' Update of bayesian parameter u
#'
#' @param A an array of dim=c(N,N,V) 
#' @param params list of parameters of the model
#'
#' @return params with u updated
#' @export
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

#' Update of bayesian parameter beta
#'
#' @param params list of parameters of the model
#'
#' @return params with beta updated
#' @export
#'
update_beta_bayesian<-function(params){
  
  tau= params$tau
  beta_0 = params$beta_0
  
  beta_k <- beta_0 + apply(tau,2,sum)
  params$beta_k <- beta_k
  return(params)
}

#' Update of bayesian parameter theta
#'
#' @param params list of parameters of the model
#'
#' @return params with theta updated
#' @export
#'
update_theta_bayesian<-function(params){
  
  u= params$u
  theta_0 = params$theta_0
  
  theta <- theta_0 + apply(u,2,sum)
  params$theta <- theta
  return(params)
}

#' Update of bayesian parameter eta
#'
#' @param A an array of dim=c(N,N,V) 
#' @param params list of parameters of the model
#'
#' @return params with eta updated
#' @export
#'
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

#' Update of bayesian parameter xi
#'
#' @param A an array of dim=c(N,N,V) 
#' @param params list of parameters of the model
#'
#' @return params with xi updated
#' @export
#'
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
