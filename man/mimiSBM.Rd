% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesianMixtureSBM.R
\name{mimiSBM}
\alias{mimiSBM}
\title{Mixture of Multilayer Integrator SBM (mimiSBM)}
\usage{
mimiSBM(
  A,
  Kset,
  Qset,
  beta_0 = 1/2,
  theta_0 = 1/2,
  eta_0 = 1/2,
  xi_0 = 1/2,
  criterion = "ILVB",
  tol = 0.001,
  iter_max = 10,
  n_init = 1,
  alternate = F,
  Verbose = F,
  eps_conv = 1e-04,
  type_init = "SBM"
)
}
\arguments{
\item{A}{an array of dim=c(N,N,V)}

\item{Kset}{Set of number of clusters}

\item{Qset}{Set of number of components}

\item{beta_0}{hyperparameters for beta}

\item{theta_0}{hyperparameters for theta}

\item{eta_0}{hyperparameters for eta}

\item{xi_0}{hyperparameters for xi}

\item{criterion}{model selection criterion, criterion=c("ILVB","ICL_approx","ICL_variationnel","ICL_exact")}

\item{tol}{convergence parameter on ELBO}

\item{iter_max}{maximal number of iteration of mimiSBM}

\item{n_init}{number of initialization of the mimi algorithm.}

\item{alternate}{boolean indicated if we put an M-step after each part of the E-step, after u optimization and after tau optimization. If not, we optimize u and tau and after the M-step is made.}

\item{Verbose}{boolean for information on model fitting}

\item{eps_conv}{parameter of convergence for tau.}

\item{type_init}{select the type of initialization type_init=c("SBM","Kmeans","random")}
}
\value{
The best model, conditionnally to the criterion, and its parameters.
}
\description{
mimiSBM is a model that allows both clustering of individuals and grouping of views by component.
This Bayesian model estimates the probability of individuals belonging to each cluster (cluster crossing all views) and the membership component for all views.
In addition, the connectivity tensor between classes, conditional on the components, is also estimated.
}
\examples{
set.seed(42)
K = c(2,3); pi_k = rep(1/4,4) ; rho = rep(1/2,2)
res <- rSMB_partition(N = 50,V = 5,K = K ,pi_k = pi_k ,rho = rho,p_switch = 0.1)
A = res$simulation$A ; Kset = 4 ; Qset = 2
model <- mimiSBM(A,Kset,Qset,n_init = 1, Verbose=FALSE)
}
