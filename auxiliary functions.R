library(coda)
library(ggplot2)
library(quantreg)
library(nimble)
library(matrixStats)

#########################
##### Distribution ######
#########################

## Generalized Gamma distribution
dGG <- function(t, scale, shape, tau, log = FALSE){
  val <- log(tau) - (shape*tau)*log(scale) - lgamma(shape) + ((shape*tau) - 1)*log(t) -
    (t/scale)^tau
  if(log) return(val) else return(exp(val))
}

pGG <- function(t, scale, shape, tau, log.p = FALSE){
  val <- pgamma( (t/scale)^tau, shape = shape, log.p = TRUE) 
  if(log.p) return(val) else return(exp(val))
}

rGG <- function(n, scale, shape, tau){
  p <- runif(n)
  out <- qgamma(p, shape = shape)^(1/tau)*scale
  return(as.vector(out))
}

## FmixUGG distribution under assumption of
## longitudinal data

dUGG <- function(x, kappa=0.2, gama=2, tau=2, q=0.5, log = FALSE) {
  n_arg <- max(length(x), length(kappa), length(gama), length(tau), length(q))
  x <- rep(x, length = n_arg)
  kappa <- rep(kappa, length = n_arg)
  gama <- rep(gama, length = n_arg)
  tau <- rep(tau, length = n_arg)
  q <- rep(q, length = n_arg)
  h<- ((1-kappa)/kappa)*(1/qgamma(1-q, gama))^(1/tau)
  d <- (1/x^2)*dGG((1-x)/x,scale=h,shape=gama,tau=tau) 
  
  if (log== F) {
    return(d)
  }
  else {
    return(log(d))
  }
}

rmixUGGl <- function(n, t, 
                     kappa=matrix(c(0.2,0.2,0.2),ncol=3), 
                     gama=c(2,2,2), 
                     theta=c(0.5,0.2,0.3), 
                     tau=c(2,2,2), 
                     q=0.5) {
  k=ncol(kappa)
  kappal=thetal=gamal=taul=h=matrix(0,ncol=k,nrow=n*t)
  for(j in 1:k){
    kappal[,j]=rep(kappa[1:nrow(kappa),j],length = n*t)
    thetal[,j]=rep(theta[j], length = n*t)
    gamal[,j]=rep(gama[j], length = n*t)
    taul[,j]=rep(tau[j], length = n*t)
    h[,j]<-(((1-kappal[,j])/kappal[,j])*
              (1/qgamma(1-q, gamal[,j]))^(1/taul[,j]))
  }
  q <- rep(q, length = n*t)
  u <- runif(n)
  out=0
  count=0
  for(i in 1:n){
    aux=cumsum(thetal[i,])
    for(j in 1:k){
      if(u[i]<aux[j]){
        index = (i+(i-1)*(t-1)):(i+i*(t-1))
        v <- rGG(t, h[index,j],gamal[index,j],taul[index,j])
        out[index]<-(1/(1+v))
        count[index]=j
        break
      }else{next}
    }
  }
  return(cbind(out,count))
}

dmixUGGl <- function(x,kappa=matrix(c(0.2,0.2,0.2),ncol=3),
                     delta=c(2,2,2), theta=c(0.5,0.2,0.3), 
                     tau=c(2,2,2), q=0.5, log = FALSE) {
  k=ncol(kappa)
  n=length(x)
  den=kappal=thetal=deltal=taul=h=matrix(0,ncol=k,nrow=n)
  for(j in 1:k){
    kappal[,j]=rep(kappa[1:nrow(kappa),j],length = n)
    thetal[,j]=rep(theta[j], length = n)
    deltal[,j]=rep(delta[j], length = n)
    taul[,j]=rep(tau[j], length = n)
    h[,j]<-((1-kappal[,j])/kappal[,j])*(1/qgamma(1-q, deltal[,j]))^(1/taul[,j])
    den[,j]<-thetal[,j] *(1/x^2)*dGG((1-x)/x,scale=h[,j],shape=deltal[,j],tau=taul[,j])
  }
  if (log== F) {
    return(rowSums(den))
  }
  else {
    return(log(rowSums(den)))
  }
}

#########################
######## Nimble #########
#########################

dggamma <- nimbleFunction(
  run = function(x = double(), scale = double(), 
                 shape = double(), tau = double(), 
                 log = logical(0, default = 0)) {
    returnType(double())
    d <- log(tau) - (shape*tau)*log(scale) - 
      lgamma(shape) + ((shape*tau) - 1)*log(x) -
      (x/scale)^tau
    if (log) return(d)
    else return(exp(d))
  })


registerDistributions(list(
  dggamma = list(
    BUGSdist = "dggamma(scale, shape, tau)",
    types = c('value = double()', 'scale = double()',
              'shape = double()', 'tau = double()')
  )))

# Half-Cauchy
dhc <- nimbleFunction(
  run = function(x = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    if(x<0){
      if(log){
        return(-Inf)
      }else{
        return(0)
      }
    }
    logProb <- log(2/(pi*(1+x^2)))
    if(log) return(logProb)
    else return(exp(logProb))
  })

rhc <- nimbleFunction(
  run = function(n = integer(0)) {
    returnType(double(0))
    if(n != 1) print("rhc only allows n = 1; using n = 1.")
    dev <- rnorm(1)/rnorm(1)
    return(abs(dev))
  })

registerDistributions(list(
  dhc = list(
    BUGSdist = "dhc()",
    range = c(0, Inf)
  )
))



#########################
###### Parameters #######
#########################

extract_param <- function(samples_matrix, 
                          param_name, 
                          dims) {
  # Filtrar colunas que correspondem ao parâmetro
  param_cols <- grep(paste0("^", param_name, "\\["), 
                     colnames(samples_matrix))
  param_data <- samples_matrix[, param_cols, 
                               drop = FALSE]
  
  # Extrair os índices dos parâmetros: de "param[2,3]" → c(2,3)
  index_matrix <- do.call(rbind, lapply(
    colnames(param_data),
    function(n) as.integer(
      unlist(regmatches(n, gregexpr("\\d+", n))))
  ))
  
  n_iter <- nrow(param_data)
  full_dims <- c(n_iter, dims)
  param_array <- array(NA, dim = full_dims)  # [n_iter x dim1 x dim2 ...]
  
  for (k in seq_len(ncol(param_data))) {
    idx <- index_matrix[k, ]
    param_array[matrix(c(1:n_iter, rep(idx, each = n_iter)), ncol = length(full_dims))] <- param_data[, k]
  }
  
  return(param_array)
}

colModes <- function(x) {
  apply(x, 2, function(col) {
    ux <- unique(col)
    ux[which.max(tabulate(match(col, ux)))]
  })
}

extract_all_params <- function(mod) {
  
  p=ncol(mod$X)
  k=mod$k
  n=mod$n
  
  samples_mat <- as.matrix(mod$samples)
  
  param_array <- extract_param(
    samples_mat, "param", c(p, k)
  )
  
  param_ale_array <- extract_param(
    samples_mat, "param_ale", c(n)
  )
  
  sigma_array <- extract_param(
    samples_mat, "sigma", c(k)
  )
  
  gama_array <- extract_param(
    samples_mat, "gama", c(k)
  )
  
  tau_array <- extract_param(
    samples_mat, "tau", c(k)
  )
  
  thetas_array <- extract_param(
    samples_mat, "thetas", c(k)
  )
  
  z_array <- extract_param(
    as.matrix(mod$samples), "z", c(n)
  )
  
  class = colModes(z_array)
  
  list(
    param_array     = param_array,
    param_ale_array = param_ale_array,
    sigma_array     = sigma_array,
    gama_array      = gama_array,
    tau_array       = tau_array,
    thetas_array    = thetas_array,
    z_array         = z_array,
    class           = class
  )
}


summarize_param <- function(samples_matrix, param_name = "beta", thetas_array = NULL) {
  dims <- dim(samples_matrix)
  
  # Detecta se tem 2 ou 3 dimensões
  if (length(dims) == 3) {
    n_samples <- dims[1]
    n_params <- dims[2]
    n_components <- dims[3]
  } else if (length(dims) == 2) {
    n_samples <- dims[1]
    n_params <- 1
    n_components <- dims[2]
    # Converte 2D para 3D: [samples, params, components]
    samples_matrix <- array(samples_matrix, dim = c(n_samples, n_params, n_components))
  } else {
    stop("Array must be 2D or 3D")
  }
  
  df_list <- list()
  
  for (k in 1:n_components) {
    # Filtra componentes pelo theta se thetas_array existir
    if (!is.null(thetas_array)) {
      theta_mean <- mean(thetas_array[, k])
      if (theta_mean <= 0.1) next
    }
    
    for (j in 1:n_params) {
      param_samples <- samples_matrix[, j, k]  # samples first
      HPD <- HPDinterval(as.mcmc(param_samples))
      
      name <- switch(param_name,
                     beta  = paste0("$\\beta_{", j-1, "}^{(", k, ")}$"),
                     sigma = paste0("$\\sigma^{(", k, ")}$"),
                     gamma = paste0("$\\gamma^{(", k, ")}$"),
                     tau   = paste0("$\\tau^{(", k, ")}$"),
                     theta = paste0("$\\theta^{(", k, ")}$"),
                     paste0(param_name, "_", j))
      
      df_list[[length(df_list) + 1]] <- data.frame(
        Parameter  = name,
        Estimate     = median(param_samples),
        "Std Deviation" = sd(param_samples),
        HPD        = paste0("[", round(HPD[1],3), ", ", round(HPD[2],3), "]")
      )
    }
  }
  
  df <- do.call(rbind, df_list)
  return(df)
}

##############################
####### Classification #######
##############################

extract_classes_filtered <- function(mod, threshold = 0.1){
  params <- extract_all_params(mod)
  
  z_array <- params$z_array
  class_vec <- params$class
  thetas_array <- params$thetas_array
  theta_means <- colMeans(thetas_array)
  valid_comps <- which(theta_means > threshold)
  class_vec[!class_vec %in% valid_comps] <- NA
  
  list(class = class_vec, thetas = theta_means,
       valid_comps = valid_comps)
}

# relabel_by_theta <- function(class_list){
#   # extrai classes e thetas
#   class_vec <- class_list$class
#   theta_means <- class_list$thetas
#   
#   # identifica apenas componentes "válidas"
#   valid_comps <- which(theta_means > 0.1)
#   
#   # ordena as componentes válidas pelo valor de theta
#   ordem <- order(theta_means[valid_comps], decreasing=TRUE)  # crescente; use decreasing=TRUE se preferir
#   nova_ordem <- valid_comps[ordem]
#   
#   # cria um mapa: componente antiga -> nova etiqueta sequencial (1, 2, 3, ...)
#   mapa <- setNames(seq_along(nova_ordem), nova_ordem)
#   
#   # aplica o mapa ao vetor de classes
#   class_vec_rerot <- ifelse(class_vec %in% valid_comps, mapa[as.character(class_vec)], NA)
#   
#   list(class = class_vec_rerot, thetas = theta_means[nova_ordem])
# }


#########################
####### Residuals #######
#########################

simulate_dharma_posterior_mean <- function(mod,
    n_sims = 1000        # número de simulações a realizar
) {
  
  y_obs=mod$y
  X=mod$X
  idv=mod$idv
  p=ncol(mod$X)
  q=mod$q
  k=mod$k
  samples_matrix=as.matrix(mod$samples)
  
  
  # Extrair médias a posteriori dos parâmetros
  param_mean     <- apply(
    extract_param(samples_matrix, "param", 
                  c(p, k)), c(2,3), mean)
  param_ale_mean <- colMeans(
    extract_param(samples_matrix, "param_ale", c(n)))
  gama_mean      <- colMeans(
    extract_param(samples_matrix, "gama", c(k)))
  thetas_mean    <- colMeans(
    extract_param(samples_matrix, "thetas", c(k)))
  tau_mean       <- colMeans(
    extract_param(samples_matrix, "tau", c(k)))
  sigma_mean     <- colMeans(
    extract_param(samples_matrix, "sigma", c(k)))
  z_mean         <- colModes(
    extract_param(samples_matrix, "z", c(n)))
  
  n_obs <- length(y_obs)
  y_sim_matrix <- matrix(NA, nrow = n_obs, ncol = n_sims)
  t <- length(unique(idv))
  
  for (s in 1:n_sims) {
    for (i in 1:n) {
      idx_i <- which(idv == i)
      t_i <- length(idx_i)
      kappa_ik <- matrix(NA, nrow = t_i, ncol = k)
      
      for (l in 1:k) {
        beta_l <- param_mean[, l]
        if (l == z_mean[i]) {
          b_il <- param_ale_mean[i]
        } else {
          b_il <- rnorm(1, 0, sqrt(sigma_mean[l]))
        }
        
        for (j in 1:t_i) {
          row_idx <- idx_i[j]
          eta_ijl <- X[row_idx, ] %*% beta_l + b_il
          kappa_ik[j, l] <- ilogit(eta_ijl)
        }
      }
      
      y_sim_matrix[idx_i, s] <- rmixUGGl(
        n = 1,
        t = t_i,
        kappa = kappa_ik,
        gama = gama_mean,
        theta = thetas_mean,
        tau = tau_mean,
        q = q
      )[, 1]
    }
  }

  scaled_residuals <- rowMeans(
    y_sim_matrix < matrix(y_obs, 
                          nrow = n_obs, 
                          ncol = n_sims),
    na.rm = TRUE
  )
  centered_residuals <- scaled_residuals - 0.5
  
  return(list(
    scaled_residuals = scaled_residuals,
    centered_residuals = centered_residuals,
    y_sim_matrix = y_sim_matrix
  ))
}


plot_qq_simple_dharma <- function(scaled_residuals) {
  
  scaled_residuals <- na.omit(scaled_residuals)
  n_obs <- length(scaled_residuals)
  # Quantis teóricos da U(0,1)
  qq_vals <- ppoints(n_obs)
  
  # Data frame para ggplot
  qq_df <- data.frame(
    theoretical = qq_vals,
    residuals = sort(scaled_residuals)
  )
  
  # QQ-plot simples
  ggplot(qq_df, aes(x = theoretical, y = residuals)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      #title = "QQ-Plot dos Resíduos Escalados",
      x = "Theoretical Quantiles U(0,1)",
      y = "Scaled Residuals"
    ) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "grey"),
      text=element_text(size=25,family="serif")
    )
}

qq_ale_plot <- function(mod, threshold, conf) {
  
  library(ggplot2)
  library(ggpubr)
  
  params <- extract_all_params(mod)
  
  param_ale_array <- params$param_ale_array
  sigma_array     <- params$sigma_array
  thetas_array    <- params$thetas_array
  class_vec       <- params$class
  
  theta_means <- colMeans(thetas_array)
  valid_comps <- which(theta_means > threshold)
  
  if (length(valid_comps) == 0) {
    warning("Nenhuma iteração satisfaz theta > threshold")
    return(NULL)
  }
  
  param_ale_array_f <- param_ale_array[valid_comps, , drop = FALSE]
  sigma_array_f     <- sigma_array[valid_comps, , drop = FALSE]
  
  qq_single <- function(ef) {
    
    cols_ef <- which(class_vec == ef)
    
    if (length(cols_ef) == 0)
      return(NULL)
    
    norm <- colMeans(param_ale_array_f[, cols_ef, drop = FALSE])
    
    if (all(is.na(norm)) || length(norm) == 0)
      return(NULL)
    
    sd_custom <- sqrt(colMeans(sigma_array_f)[ef])
    
    if (is.na(sd_custom) || sd_custom == 0)
      return(NULL)
    
    ef_ale <- data.frame(norm_std = norm / sd_custom)
    
    ggplot(ef_ale, aes(sample = norm_std)) +
      stat_qq_band(conf = conf) +
      stat_qq_line() +
      stat_qq_point() +
      labs(
        title = paste("Componente", ef),
        x = "Theoretical Quantiles",
        y = "Quantile Residuals"
      ) +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey"),
        text = element_text(size = 20, family = "serif")
      )
  }
  
  plots <- lapply(seq_len(ncol(sigma_array_f)), qq_single)
  plots <- Filter(Negate(is.null), plots)
  
  if (length(plots) == 0) {
    warning("Nenhum gráfico válido para plotar")
    return(NULL)
  }

  ggarrange(
    plotlist = plots,
    ncol = 2,
    nrow = ceiling(length(plots) / 2)
  )
}


################################
#### Information Criteria ######
################################

calc_criteria <- function(mod) {
  
  X <- mod$X
  y <- mod$y
  idv <- mod$idv
  n <- mod$n
  t <- mod$t
  k <- mod$k
  q <- mod$q
  
  params <- extract_all_params(mod)
  param_array     <- params$param_array
  param_ale_array <- params$param_ale_array
  sigma_array     <- params$sigma_array
  gama_array      <- params$gama_array
  tau_array       <- params$tau_array
  z_array         <- params$z_array
  class_vec       <- params$class
  
  n_obs <- n * t
  m <- dim(param_array)[1]
  
  # ======================
  #        WAIC
  # ======================
  lppd <- pwaic <- numeric(n_obs)
  
  for (i in 1:n_obs) {
    ind_i <- idv[i]
    den <- numeric(m)
    
    for (j in 1:m) {
      z_ji   <- z_array[j, ind_i]
      beta_l <- param_array[j, , z_ji]
      b_il   <- param_ale_array[j, ind_i]
      
      xb <- sum(X[i, ] * beta_l) + b_il
      kappa <- max(min(ilogit(xb), 1 - 1e-4), 1e-4)
      
      den[j] <- dUGG(x = y[i],
                     kappa = kappa,
                     gama  = gama_array[j, z_ji],
                     tau   = tau_array[j, z_ji],
                     q     = q,
                     log   = FALSE)
    }
    
    lppd[i]  <- log(mean(den))
    pwaic[i] <- var(log(den))
  }
  
  WAIC <- -2 * (sum(lppd) - sum(pwaic))
  
  # ======================
  #        LPML
  # ======================
  aux <- matrix(0, m, n)
  
  for (j in 1:m) {
    for (i in 1:n) {
      z_i <- z_array[j, i]
      dens_i <- 1
      
      for (tt in 1:t) {
        obs_index <- (i - 1) * t + tt
        beta_l <- param_array[j, , z_i]
        b_il   <- param_ale_array[j, i]
        
        xb <- sum(X[obs_index, ] * beta_l) + b_il
        kappa <- max(min(ilogit(xb), 1 - 1e-4), 1e-4)
        
        dens_tt <- dUGG(x = y[obs_index],
                        kappa = kappa,
                        gama  = gama_array[j, z_i],
                        tau   = tau_array[j, z_i],
                        q     = q,
                        log   = FALSE)
        
        dens_i <- dens_i * dens_tt
      }
      
      aux[j, i] <- dens_i
    }
  }
  
  cpo  <- 1 / colMeans(1 / aux)
  LPML <- sum(log(cpo))
  
  # ======================
  #         DIC
  # ======================
  
  ale_hat    <- colMedians(param_ale_array)
  z_mean    <- colModes(extract_param(as.matrix(mod$samples), "z", c(n)))
  
  lvero_vec <- numeric(n_obs)
  
  for (i in 1:n_obs) {
    id_i <- idv[i]
    z_i  <- z_mean[id_i]
    b_il <- ale_hat[id_i]
    lobs <- 0
    
    beta_l <- colMedians(param_array[, , z_i])
    gama_i <- median(gama_array[, z_i])
    tau_i  <- median(tau_array[, z_i])
    
    for (tt in 1:t) {
      obs_index <- (id_i - 1) * t + tt
      
      xb <- sum(X[obs_index, ] * beta_l) + b_il
      kappa <- max(min(ilogit(xb), 1 - 1e-4), 1e-4)
      
      lobs <- lobs + dUGG(x = y[obs_index],
                          kappa = kappa,
                          gama  = gama_i,
                          tau   = tau_i,
                          q     = q,
                          log   = TRUE)
    }
    
    lvero_vec[i] <- lobs
  }
  
  lvero  <- sum(lvero_vec)
  BARDEV <- mean(colSums(-2 * log(aux)))
  DEVBAR <- -2 * lvero
  DIC    <- DEVBAR + 2 * (BARDEV - DEVBAR)
  
  # ======================
  #        OUTPUT
  # ======================
  list(
    WAIC = WAIC,
    LPML = LPML,
    DIC  = DIC
  )
}



