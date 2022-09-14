## ----install pacakges, eval=FALSE-----------------------------------------------
## # For Windows:
## install.packages("survHE",
## 	               repos=c("http://www.statistica.it/gianluca/R",
## 		                     "https://cran.rstudio.org",
##                          "https://inla.r-inla-download.org/R/stable"),
## 	               dependencies=TRUE)
## 
## install.packages("INLA",
##                  repos=c(getOption("repos"),
##                          INLA="https://inla.r-inla-download.org/R/stable"),
##                  dep=TRUE)


## ----load packages, results='hide', warning=FALSE, message=FALSE----------------
library(survHE)
library(INLA)


## ----load data------------------------------------------------------------------
load("../Data/ta174_FCR.Rdata")
head(dat_FCR)


## ----source function, include=TRUE----------------------------------------------
source("../Scripts/Functions.R", local = knitr::knit_global())


## ----observed estimate----------------------------------------------------------
obs_Surv <- surv_est_inla(data = dat_FCR, 
                          cutpoints = seq(0,180, by =5))
                          
names(obs_Surv)


## ----observedplot, echo=TRUE, fig.cap="Survival estimates of the FCR arm based on the piecewise exponential model fitted to the observed data"----
ggplot() + geom_line(aes(obs_Surv$KM$time, obs_Surv$KM$surv, colour = "Kaplan-Meier"), size = 1.25, linetype = "dashed")+
  xlim(0, 180) + ylim(0,1) +
  geom_line(aes(obs_Surv$time,rowMeans(obs_Surv$S_obs), colour = "Data fitting"), size = 1)+
  geom_ribbon(aes(x=obs_Surv$time, y = rowMeans(obs_Surv$S_obs),
                  ymin = apply(obs_Surv$S_obs, 1, quantile, probs = 0.025),ymax = apply(obs_Surv$S_obs, 1, quantile, probs = 0.975)), alpha = 0.1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), text = element_text(size = 8))+
  xlab("Time (months)") + ylab("Survival")+
  scale_colour_manual(name = "model", values = c("Data fitting"="#7CAE00", "Kaplan-Meier"="brown"))+
  annotate("text", x = 150, y = 0.9, label = "FCR arm", colour = "orange", fontface = 2)


## ----external estimate, message=FALSE, warning=FALSE, results='hide'------------
ext_Surv <- ext_surv_est(t_info = 144, S_info = 0.05, 
                         T_max = 180, 
                         times_est = seq(0, 180), distr = "gom")

names(ext_Surv)


## ----external_object, echo=FALSE------------------------------------------------
names(ext_Surv)


## ----weight function------------------------------------------------------------
blending_interval <- list(min = 48, max = 150)
beta_params <- list(alpha = 3, beta = 3)

tp <- seq(0, 180)

# parameters for the weight function
wt_par <- list(a = blending_interval$min, b = blending_interval$max,
               shape1 = beta_params$alpha, shape2 = beta_params$beta) 

weight <- with(wt_par,
               pbeta((tp - a)/(b - a), shape1, shape2))


## ----blended estimate-----------------------------------------------------------
n_sim <- ncol(ext_Surv$S_ext)

ble_Surv <- matrix(NA, nrow = 180 + 1, ncol = n_sim)

for (i in seq_len(n_sim)) {
  ble_Surv[, i] <-
    obs_Surv$S_obs[, i]^(1 - weight) * ext_Surv$S_ext[, i]^weight
}



## ----blendedplot, echo=TRUE, fig.cap="Blended survival curve based on short-term data and external information for FCR arm"----
ggplot() + geom_line(aes(obs_Surv$KM$time, obs_Surv$KM$surv, colour = "Kaplan-Meier"), size = 1.25, linetype = "dashed")+
  xlim(0, 180) + ylim(0,1) +
  geom_line(aes(tp,rowMeans(obs_Surv$S_obs), colour = "Data fitting"), size = 1, linetype = "twodash")+
  geom_ribbon(aes(x=tp, y = rowMeans(obs_Surv$S_obs),
                  ymin = apply(obs_Surv$S_obs, 1, quantile, probs = 0.025),
                  ymax = apply(obs_Surv$S_obs, 1, quantile, probs = 0.975)), alpha = 0.1)+
  geom_line(aes(tp, rowMeans(ext_Surv$S_ext), colour = "External info"), size = 1, linetype = "longdash")+
  geom_ribbon(aes(x=tp, y = rowMeans(ext_Surv$S_ext), ymin = apply(ext_Surv$S_ext, 1, quantile, probs = 0.025),
                  ymax = apply(ext_Surv$S_ext, 1, quantile, probs = 0.975)), alpha = 0.1)+
  geom_line(aes(tp, rowMeans(ble_Surv), colour = "Blended curve"), size = 1.25)+
  geom_ribbon(aes(x=tp, y = rowMeans(ble_Surv), ymin = apply(ble_Surv, 1, quantile, probs = 0.025),
                  ymax = apply(ble_Surv, 1, quantile, probs = 0.975)), alpha = 0.1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black"), text = element_text(size = 8))+
  xlab("Time (months)") + ylab("Survival")+
  scale_colour_manual(name = "model", values = c("Data fitting"="#7CAE00", "External info"="#00BFC4", "Blended curve"="#F8766D", "Kaplan-Meier"="brown"))+
  annotate("text", x = 150, y = 0.9, label = "FCR arm", colour = "orange", fontface = 2)


