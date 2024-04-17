library(cmdstanr)
library(posterior)
library(dqrng)

###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###             BRCAL: Bayesian Regularized Calibrated              ###
###               Posterior draws of ATE via BRCAL                  ###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###


# Input -------------------------------------------------------------------
# Y: Outcome vector
# A: Treatment vector
# X: Matrix of covariates
# weight: learning rate of general posterior
# PCIC.print: if TRUE, output the Posterior Covariance Information Criterion (default TRUE) 
# mu1, tau1: mean and precision prior parameter of potential outcome if exposured
# mu0, tau0: mean and precision prior parameter of potential outcome if exposured
# iter_warmup: length of burn-in
# iter_sampling: length of posterior draws
# thin: period between saved samples

# Output (vector) ---------------------------------------------------------
# $post_ate_brcal: posterior draws of ATE 
# $PCIC: posterior covariance information criterion
# $post_lambda: posterior draws of a specified upper bound for covariate balance


# Loading stan file -------------------------------------------------------
stan_ps_brcal <- cmdstanr::cmdstan_model("ps_brcal.stan")

# Main function -----------------------------------------------------------
BRCAL <- function(Y, A, X, weight=1, PCIC.print=TRUE,
                  mu1=0, mu0=0, tau1=100^(-2), tau0=100^(-2),
                  iter_warmup=1000, iter_sampling=2000, thin=1){
  dlist <- list(N=dim(X)[1], D=dim(X)[2], A=A, X=X, R=rep(1,dim(X)[1]), w_alpha=weight)
  
  ## sampling of propensity score -------------------------------------------
  fit_alpha <- stan_ps_brcal$sample(data=dlist, seed=1234, refresh=0, init=0,
                                    chains=1, iter_warmup=iter_warmup, iter_sampling=iter_sampling, thin=thin,
                                    show_messages=FALSE)
  post_alpha <- posterior::as_draws_matrix(fit_alpha$draws(c("alpha0","alpha")))
  post_ps_brcal <- (1+exp(-tcrossprod(post_alpha,cbind(1,X))))^{-1}
  
  ## sampling of ATE --------------------------------------------------------
  s1 <- t(apply(post_ps_brcal, 1, function(x){(2*A)/x}))
  s0 <- t(apply(post_ps_brcal, 1, function(x){(2*(1-A))/(1-x)}))
  til_tau1 <- tau1+weight*rowSums(s1)
  til_tau0 <- tau0+weight*rowSums(s0)
  til_mu1 <- til_tau1^{-1} * ( tau1*mu1 + weight*apply(s1,1,function(x){sum(x*Y)}) )
  til_mu0 <- til_tau0^{-1} * ( tau0*mu0 + weight*apply(s0,1,function(x){sum(x*Y)}) )
  
  ## calculate PCIC ---------------------------------------------------------
  if(PCIC.print) {
    post_theta1 <- til_mu1+sqrt(til_tau1^(-1))*dqrng::dqrnorm(iter_sampling)
    post_theta0 <- til_mu0+sqrt(til_tau0^(-1))*dqrng::dqrnorm(iter_sampling)
    
    post_ate_brcal <- post_theta1 - post_theta0

    n <- length(Y)
    v_theta1 <- -t(apply(post_ps_brcal,1,function(x){A/x}))*(t(matrix(rep(Y,iter_sampling),n))-post_theta1)^2
    v_theta0 <- -t(apply(post_ps_brcal,1,function(x){(1-A)/(1-x)}))*(t(matrix(rep(Y,iter_sampling),n))-post_theta0)^2
    v_alpha <- -A*exp(-tcrossprod(post_alpha,cbind(1,X))) - (1-A)*tcrossprod(post_alpha,cbind(1,X)) - (1-A)*exp(tcrossprod(post_alpha,cbind(1,X))) + A*tcrossprod(post_alpha,cbind(1,X))
    
    TG <- mean({v_theta1+v_theta0+v_alpha})
    V <- mean(({v_theta1+v_theta0+v_alpha})^2) - mean((rowMeans({v_theta1+v_theta0+v_alpha}))^2)
    PCIC <- TG - V
    
    return(list(post_ate_brcal=post_ate_brcal, post_lambda=posterior::as_draws_matrix(sqrt(fit_alpha$draws("lambda"))), PCIC=PCIC))
  }else {
    sig <- sqrt( (til_tau1)^{-1} + (til_tau0)^{-1} )
    mu <- til_mu1 - til_mu0
    post_ate_brcal <- mu + sig*dqrng::dqrnorm(iter_sampling)
    return(post_ate_brcal)
  }
  
}
