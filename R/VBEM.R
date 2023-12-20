#---------------- Variational Bayes EM ----------------

#' Variational Bayes Expectation Maximization
#'
#' @param A an array of dim=c(N,N,V) 
#' @param params list of parameters of the model
#' @param alternate boolean indicated if we put an M-step after each part of the E-step, after u optimization and after tau optimization. If not, we optimize u and tau and after the M-step is made.
#' @param eps_conv parameter of convergence for tau.
#'
#' @return  params with updated parameters.
#' @export
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