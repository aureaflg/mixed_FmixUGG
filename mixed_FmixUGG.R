library(coda)
library(ggplot2)
library(qqplotr)
library(ggpubr)
library(quantreg)
library(nimble)
library(matrixStats)
library(knitr)

source("auxiliary functions.R")


###################### Function to fit the Mixed FmixUGG regression model #######################

## With nchains you set the number of MCMC chains, with thin the user can
## specify how much the MCMC chains should be thinned out before storing them. niter
## is used to set the iteraction numbers and nburnin is used to set the burn-in.
## The argument q defines which quantile will be modelled (0<q<1).  

FmixUGG=function(kappa.formula=formula,
                 data=NULL, time=NULL, k=3,
                 q=0.5, nchains = 1,thin=50, 
                 niter = 90000, nburnin = 40000){
  
  mcall <- if(is.null(data)){ terms(kappa.formula) 
  }else terms(kappa.formula, data = data)
  X <- if(is.null(data)){ model.matrix(kappa.formula) 
  }else model.matrix(kappa.formula, data = data)  
  y <- if(is.null(data)){ model.frame(mcall)[,1] 
  }else model.frame(mcall, data = data)[,1]
  
  if(is.null(data)){
    t=time
  }else{
    t=length(unique(data$time))
  }
  n=length(y)/t
  p=ncol(X)
  idv=rep(1:n, each=t)
  
  Code <- nimbleCode({
    for(j in 1:k){
      for(s in 1:p){
        param[s,j] ~ dnorm(0, tau=1/1000)
      }
      psi1[j] ~ dunif(0.08,1)
      psi2[j] ~ dunif(0.08,1)
      sigma[j] ~ dhc()
      gama[j] ~ dexp(psi1[j])
      tau[j] ~ dexp(psi2[j])
    }
    alpha ~ dgamma(1, 1)
    for(b in 1:(k-1)){
      v[b] ~ dbeta(1, alpha)
    }
    thetas[1:k] <- stick_breaking(v[1:(k-1)])
    for(i in 1:n){
      z[i] ~ dcat(thetas[1:k])
      param_ale[i] ~ dnorm(0, tau = 1 / sigma[z[i]])
    }
    for(i in 1:ntotal){
      for(j in 1:k){
        a[i,j] <- inprod(X[i,1:p],param[1:p,j])+
          param_ale[idv[i]]
        kappa[i,j] <- ilogit(a[i,j])
        h[i,j]<- (((1-kappa[i,j])/kappa[i,j])*
                    (1/qgamma(1-q, gama[j]))^(1/tau[j]))
        l[i,j] <- (log((1/y[i]^2))+dggamma((1-y[i])/y[i], 
                                           shape = gama[j], 
                                           scale=h[i,j],
                                           tau=tau[j], 
                                           log=TRUE))
      }
    }
    for(o in 1:n){
      density[o] <- C - sum(l[((o-1)*t+1):(o*t), z[o]])
      zeros[o] ~ dpois(density[o])
    }
  })
  
  Consts <- list(n = n, p = p, t = t, 
                 ntotal = n*t, q = q, k = k, 
                 C=10000, 
                 idv = idv)
  
  Inits<- list(param = matrix(0.1,ncol=k,nrow=p), 
               thetas = rdirch(1,rep(1/k,k)), 
               sigma = rep(1,k), z = rep(1,n),
               gama = rep(1,k),
               tau = rep(1,k),
               psi1 = rep(0.1,k),
               psi2 = rep(0.1,k),
               param_ale = rep(0.01,n),
               v=rbeta(k, 1, 1),
               alpha=1)
  
  
  Data <- list(y = y, X=X, zeros=rep(0,n))
  
  mcmc.out <- nimbleMCMC(code = Code, constants = Consts,
                         data = Data, inits = Inits,
                         nchains = nchains,thin=thin, 
                         niter = niter, 
                         nburnin = nburnin,
                         monitors = c('param',
                                      'param_ale',
                                      'tau',
                                      'gama',
                                      'thetas',
                                      'sigma',
                                      'z'),
                         summary = TRUE)
  
  return(list(mcmc.out=mcmc.out,X=X,k=k,
              y=y,q=q,n=n,t=t,idv=idv))
}

modelsummary=function(mod){
  params <- extract_all_params(mod)
  
  param_array     <- params$param_array
  param_ale_array <- params$param_ale_array
  sigma_array     <- params$sigma_array
  gama_array      <- params$gama_array
  tau_array       <- params$tau_array
  thetas_array    <- params$thetas_array
  
  df_beta <- summarize_param(param_array, "beta", thetas_array)
  df_sigma <- summarize_param(sigma_array, "sigma", thetas_array)
  df_gamma <- summarize_param(gama_array, "gamma", thetas_array)
  df_tau <- summarize_param(tau_array, "tau", thetas_array)
  df_theta <- summarize_param(thetas_array, "theta", thetas_array)
  
  df_all <- rbind(df_beta, df_sigma, df_gamma, df_tau, df_theta)
  
  kable(
    df_all,
    format  = "pipe",
    caption = "Posterior summaries (median, variance, 95% HPD) for parameters of components with θ > 0.1"
  )
  
}

# This function returns the proportion of observations classified in each component
# as well as the classification of each observation according to the fitted model.
# Use the "threshold argument" to fix a cutoff value below which mixture components 
# with low estimated mixing weights are excluded from the model classification analysis.

classification <- function(mod, threshold = 0.1){
  
  class_list <- extract_classes_filtered(mod, threshold = threshold)
  
  yclass      <- class_list$class
  theta_means <- class_list$thetas
  valid_comps <- class_list$valid_comps
  
  k_all <- seq_along(theta_means)
  
  prop_obs_all <- prop.table(
    table(factor(yclass, levels = k_all))
  )
  
  prop_obs  <- prop_obs_all[valid_comps]
  theta_val <- theta_means[valid_comps]
  
  proportion <- data.frame(
    Component = paste("Component", valid_comps),
    "Proportion of Observations Assigned to Component" =
      as.numeric(prop_obs),
    "Estimated mixing weight" =
      theta_val,
    check.names = FALSE
  )
  
  # class_r <- relabel_by_theta(class)
  # k <- length(class_r$thetas)
  # 
  # yclass <- class_r$class
  # 
  # prop_obs <- prop.table(table(factor(yclass, levels = 1:k)))
  # 
  # proportion <- data.frame(
  #   Component = paste("Component", 1:k),
  #   "Proportion of Observations Assigned to Component" = as.numeric(prop_obs),
  #   "Estimated mixing weight" = class_r$thetas,
  #   check.names = FALSE
  # )
  
  list(
    proportion = proportion,
    "classification of observations" = yclass
  )
}

modelresiduals=function(mod){
  res = simulate_dharma_posterior_mean(mod, n_sims = 1000)
  fig1 = plot_qq_simple_dharma(res$scaled_residuals)
  df_res <- data.frame(
    Index = seq_along(res$scaled_residuals),
    Scaled_Residuals = res$scaled_residuals
  )
  
  fig2 = ggplot(df_res, aes(x = Index, y = Scaled_Residuals)) +
    geom_point(size = 2) +
    labs(
      x = "Index",
      y = "Scaled Residuals"
    ) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text = element_text(size = 25, family = "serif")
    )
  fig = ggarrange(fig1, fig2)
  print(fig)
}

# Use the "threshold argument" to fix a cutoff value below which mixture components 
# with low estimated mixing weights are excluded from the model classification analysis.

modelrandomeffects=function(mod, threshold = 0.1, conf = 0.95){
  
  fig <- qq_ale_plot(mod, threshold, conf)
  print(fig)

}

information_criteria=function(mod){
  
  ic <- calc_criteria(mod)
  
  ic_tbl <- data.frame(
    WAIC = ic$WAIC,
    LPML = ic$LPML,
    DIC  = ic$DIC
  )
  
  kable(
    ic_tbl,
    digits = 3,
    format  = "pipe",
    align = "c",
    caption = "Information Criteria"
  )
}
