#' A function to generate survival estimates with a piecewise Cox model (using INLA)
#'
#' @param inla.formula   The formula for PEM which must be a inla.surv object
#' @param dat            A dataframe for survival data with time (death_t) and event (death)
#' @param cutpoints      A sequence of cutpoints for intervals in the baseline hazard 
#' @param nsim           The number of simulations from posteriors 
#'
#' @return 
#' 


surv_est_inla <- function(inla.formula = inla.surv(death_t, death) ~ -1, data,
                          cutpoints, nsim = 100){
  
  # Convert a Cox proportional hazard model into Poisson regression
  p <- INLA::inla.coxph(
    inla.formula,
    data = data,
    control.hazard = list(
      constr = FALSE,
      cutpoints = cutpoints,
      model = "rw1"))
  
  # Fit the model 
  m <- INLA::inla(
    p$formula,
    family = p$family,
    data = c(as.list(p$data), p$data.list),
    E = p$E,
    control.compute = list(config = TRUE, dic = TRUE))
  
  # Draw samples from the joint posterior distribution  
  jp <-
    INLA::inla.posterior.sample(
      n=nsim,
      m,
      selection=list(
        Predictor=-c(1:nrow(p$data)),
        baseline.hazard=c(1:(nrow(m$summary.random$baseline.hazard))))
    )
  
  sim <-
    lapply(jp,function(i) i$latent) |>
    unlist() |>
    matrix(nrow=nsim,byrow = TRUE) |>
    as_tibble(.name_repair=~vctrs::vec_as_names(rownames(jp[[1]]$latent),quiet=TRUE))
  
  interval.t <- m$summary.random$baseline.hazard$ID
  
  # Transform baseline hazard 
  logh0 <- select(sim, contains("baseline"))
  
  if(interval.t[2] < 1){
    H0 <- (apply(exp(logh0) |>
                   select(.,1:length(interval.t)),1,cumsum))*interval.t[2]
    
    tt <- interval.t+interval.t[2]
    
  }else{
    h0 <-
      exp(logh0) |>
      select(1:length(interval.t))
    
    h0_long <- h0[,rep(1:(length(interval.t)-1), each = interval.t[2])]
    h0_long <- cbind(h0_long, h0[,length(interval.t)])
    H0= apply(h0_long,1,cumsum)
    
    tt <- 1:(max(interval.t)+1)
    
  }
  
  # Transform survival probabilities 
  S0 <- t(exp(-t(H0)))
  
  row.names(S0) <- NULL
  
  S_obs <- rbind(rep(1, nsim), S0[-nrow(S0),])
  
  # Kaplan-Meier estimate
  km <- survfit(Surv(death_t, death)~1, data = data)
  
  return(
    list(
      time = tt-tt[1],      
      S_obs = S_obs,
      KM = km))
}


#' A function to create an external survival curve based on the expert opinion 
#'
#' @param t_info   A vector of times for which expert opinion is elicited 
#' @param S_info   A vector of mean survival probabilities estimated by experts corresponding to timepoints in the t_pri  
#' @param T_max   The maximum survival time to be used 
#' @param n       The number of patients to construct the artificial external dataset 
#' @param tp      A vector of times for which the survival curves are to be computed
#' @param nsim    The number of simulations from the distribution of the survival curves
#'
#' @return
#' 


ext_surv_est <- function(t_info, S_info, T_max, times_est, 
                         n = 70, nsim = 100, distr = "gom"){ 
  
  # Sets seed for replicability
  set.seed(1996)
  
  # Partition the time horizon into intervals
  S <- c(1, S_info, 0)
  t <- c(0, t_info, T_max)
  
  c <- length(S)-1
  d <- length(t)
  
  n_par <- vector(mode = "numeric", length = c)
  
  for (i in 1:c){
    
    n_par[i] <- round(n*(S[i]-S[i+1]), digits = 0)
  }
  
  n_sim <- sum(n_par)
  
  min_unif <- rep(t[-d], n_par)
  max_unif <- rep(t[-1], n_par)
  
  # Create survival times using uniform distribution
  time   <- runif(n_sim, min_unif, max_unif)
  
  status <- rep(1, n_sim)
  
  dat_sim <- data.frame(time = time, event = status)
  
  # Fit a model to the synthetic dataset 
  if (distr == "gom"){   # Gompertz distribution 
    
    m_sim <- survHE::fit.models(formula = Surv(time, event)~1, data = dat_sim,
                                distr=distr, method="hmc",
                                priors=list(gom=list(a_alpha=0.1,b_alpha=0.1)))
    
  }else {
    
    m_sim <- survHE::fit.models(formula = Surv(time, event)~1, data = dat_sim,
                                distr = distr)
  }
  
  # Calculate the survival curve
  extr <- survHE::make.surv(m_sim,t=times_est,nsim=nsim)
  
  ext_est <- as.matrix(extr$mat[[1]])[,-1]
  
  return(list(S_ext = ext_est))
}
