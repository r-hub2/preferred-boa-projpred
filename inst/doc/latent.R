params <-
list(EVAL = TRUE)

## ----SETTINGS-knitr, include=FALSE--------------------------------------------
stopifnot(require(knitr))
knitr::opts_chunk$set(
  dev = "png",
  dpi = 150,
  fig.asp = 0.618,
  fig.width = 7,
  out.width = "90%",
  fig.align = "center",
  comment = NA,
  eval = if (isTRUE(exists("params"))) params$EVAL else FALSE
)

## ----dat_poiss----------------------------------------------------------------
# Number of observations:
N <- 100
# Data-generating function:
sim_poiss <- function(nobs = N, ncon = 4, ncats = 4, nnoise = 20) {
  # Regression coefficients for continuous predictors:
  coefs_con <- seq(from = 0.2, to = 0.8, length.out = ncon)
  # Continuous predictors:
  dat_sim <- matrix(rnorm(nobs * ncon), ncol = ncon)
  # Start linear predictor:
  linpred <- 1.4 + dat_sim %*% coefs_con
  
  # Categorical predictor:
  dat_sim <- data.frame(
    x = dat_sim,
    xcat = gl(n = ncats, k = nobs %/% ncats, length = nobs,
              labels = paste0("cat", seq_len(ncats)))
  )
  # Regression coefficients for the categorical predictor:
  coefs_cat <- seq(from = 0.2, to = 0.8, length.out = ncats)
  # Continue linear predictor:
  linpred <- linpred + coefs_cat[dat_sim$xcat]
  
  # Noise predictors:
  dat_sim <- data.frame(
    dat_sim,
    xn = matrix(rnorm(nobs * nnoise), ncol = nnoise)
  )
  
  # Poisson response, using the log link (i.e., exp() as inverse link):
  dat_sim$y <- rpois(nobs, lambda = exp(linpred))
  # Shuffle order of observations:
  dat_sim <- dat_sim[sample.int(nobs), , drop = FALSE]
  # Drop the shuffled original row names:
  rownames(dat_sim) <- NULL
  return(dat_sim)
}
# Generate data:
set.seed(300417)
dat_poiss <- sim_poiss()

## ----rstanarm_attach, message=FALSE-------------------------------------------
library(rstanarm)

## ----ref_fit_poiss------------------------------------------------------------
# Number of regression coefficients:
( D <- sum(grepl("^x", names(dat_poiss))) )
# Prior guess for the number of relevant (i.e., non-zero) regression
# coefficients:
p0 <- 5
# Prior guess for the overall magnitude of the response values, see Table 1 of
# Piironen and Vehtari (2017, DOI: 10.1214/17-EJS1337SI):
mu_prior <- 10
# Hyperprior scale for tau, the global shrinkage parameter:
tau0 <- p0 / (D - p0) / sqrt(mu_prior) / sqrt(N)
# Set this manually if desired:
ncores <- parallel::detectCores(logical = FALSE)
### Only for technical reasons in this vignette (you can omit this when running
### the code yourself):
ncores <- min(ncores, 2L)
###
options(mc.cores = ncores)
# All other columns of the data are predictors:
refm_fml <- y ~ .
refm_fit_poiss <- stan_glm(
  formula = refm_fml,
  family = poisson(),
  data = dat_poiss,
  prior = hs(global_scale = tau0, slab_df = 100, slab_scale = 1),
  ### Only for the sake of speed (not recommended in general):
  chains = 2, iter = 1000,
  ### In general, use defaults or explicitly define
  ### chains = 4, iter = 4000,
  ###
  refresh = 0
)

## ----projpred_attach, message=FALSE-------------------------------------------
library(projpred)

## ----cvvs_lat-----------------------------------------------------------------
refm_poiss <- get_refmodel(refm_fit_poiss, latent = TRUE)
time_lat <- system.time(cvvs_lat <- cv_varsel(
  refm_poiss,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0,
  ###
  seed = 95930
))

## ----time_lat-----------------------------------------------------------------
print(time_lat)

## ----plot_vsel_lat, fig.alt="Plot of the ratio of the geometric mean predictive density (GMPD) of each submodel to that of the reference model, with uncertainty bars, for the latent projection of the Poisson model. The top x-axis shows the submodel size and the bottom x-axis the predictor added at each size. A dashed horizontal line at 1 marks the reference model. The ratio increases with the size, reaches the reference model and then levels off."----
ggplot2::theme_set(bayesplot::theme_default(base_family = "sans"))
options(projpred.plot_vsel_size_position = "secondary_x")
( gg_lat <- plot(cvvs_lat, stats = "gmpd", deltas = TRUE) )

## ----size_man_lat-------------------------------------------------------------
size_decided_lat <- 5

## ----size_sgg_lat-------------------------------------------------------------
suggest_size(cvvs_lat, stat = "gmpd")

## ----predictors_final_lat-----------------------------------------------------
rk_lat <- ranking(cvvs_lat)
( predictors_final_lat <- head(rk_lat[["fulldata"]], size_decided_lat) )

## ----cvvs_trad----------------------------------------------------------------
time_trad <- system.time(cvvs_trad <- cv_varsel(
  refm_fit_poiss,
  ### Only for the sake of speed (not recommended in general):
  nclusters_pred = 20,
  nterms_max = 10,
  ###
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0,
  ###
  seed = 95930
))

## ----post_vs_trad-------------------------------------------------------------
rbind(latent = time_lat, traditional = time_trad)[, "elapsed"]
suggest_size(cvvs_trad, stat = "gmpd")
head(ranking(cvvs_trad)[["fulldata"]], size_decided_lat)

## ----plot_vsel_trad, fig.alt="Plot of the ratio of the geometric mean predictive density (GMPD) of each submodel to that of the reference model, with uncertainty bars, for the traditional projection of the Poisson model. The top x-axis shows the submodel size and the bottom x-axis the predictor added at each size. A dashed horizontal line at 1 marks the reference model. The ratio increases with the size, reaches the reference model and then levels off."----
( gg_trad <- plot(cvvs_trad, stats = "gmpd", deltas = TRUE) )

## ----ref_fit_nebin------------------------------------------------------------
refm_fit_nebin <- stan_glm(
  formula = refm_fml,
  family = neg_binomial_2(),
  data = dat_poiss,
  prior = hs(global_scale = tau0, slab_df = 100, slab_scale = 1),
  ### Only for the sake of speed (not recommended in general):
  chains = 2, iter = 1000,
  ###
  refresh = 0
)

## ----cvvs_nebin---------------------------------------------------------------
refm_prec <- as.matrix(refm_fit_nebin)[, "reciprocal_dispersion", drop = FALSE]
latent_ll_oscale_nebin <- function(ilpreds,
                                   dis = rep(NA, nrow(ilpreds)),
                                   y_oscale,
                                   wobs = rep(1, ncol(ilpreds)),
                                   cens,
                                   cl_ref,
                                   wdraws_ref = rep(1, length(cl_ref))) {
  y_oscale_mat <- matrix(y_oscale, nrow = nrow(ilpreds), ncol = ncol(ilpreds),
                         byrow = TRUE)
  wobs_mat <- matrix(wobs, nrow = nrow(ilpreds), ncol = ncol(ilpreds),
                     byrow = TRUE)
  refm_prec_agg <- cl_agg(refm_prec, cl = cl_ref, wdraws = wdraws_ref)
  ll_unw <- dnbinom(y_oscale_mat, size = refm_prec_agg, mu = ilpreds, log = TRUE)
  return(wobs_mat * ll_unw)
}
latent_ppd_oscale_nebin <- function(ilpreds_resamp,
                                    dis_resamp = rep(NA, nrow(ilpreds_resamp)),
                                    wobs = rep(1, ncol(ilpreds_resamp)),
                                    cl_ref,
                                    wdraws_ref = rep(1, length(cl_ref)),
                                    idxs_prjdraws) {
  refm_prec_agg <- cl_agg(refm_prec, cl = cl_ref, wdraws = wdraws_ref)
  refm_prec_agg_resamp <- refm_prec_agg[idxs_prjdraws, , drop = FALSE]
  ppd <- rnbinom(prod(dim(ilpreds_resamp)), size = refm_prec_agg_resamp,
                 mu = ilpreds_resamp)
  ppd <- matrix(ppd, nrow = nrow(ilpreds_resamp), ncol = ncol(ilpreds_resamp))
  return(ppd)
}
refm_nebin <- get_refmodel(refm_fit_nebin, latent = TRUE,
                           latent_ll_oscale = latent_ll_oscale_nebin,
                           latent_ppd_oscale = latent_ppd_oscale_nebin)
cvvs_nebin <- cv_varsel(
  refm_nebin,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0
  ###
)

## ----plot_vsel_nebin, fig.alt="The same kind of plot for the latent projection of the negative binomial model. The ratio increases with the size, rises slightly above 1 and then slowly decreases towards 1. The first 5 predictors of the ranking are the 5 truly relevant ones."----
( gg_nebin <- plot(cvvs_nebin, stats = "gmpd", deltas = TRUE) )

## ----size_man_nebin-----------------------------------------------------------
size_decided_nebin <- 5

## ----size_sgg_nebin-----------------------------------------------------------
suggest_size(cvvs_nebin, stat = "gmpd")

## ----predictors_final_nebin---------------------------------------------------
rk_nebin <- ranking(cvvs_nebin)
( predictors_final_nebin <- head(rk_nebin[["fulldata"]],
                                 size_decided_nebin) )

## ----data_surv----------------------------------------------------------------
N_surv <- 500
n_pred <- 50
n_pred_truth <- 10
n_pred_noise <- n_pred - n_pred_truth
dat_sim_surv <- matrix(rnorm(N_surv * n_pred), ncol = n_pred)
colnames(dat_sim_surv) <- paste0("x", c(paste0(".", seq_len(n_pred_truth)),
                                        paste0("n.", seq_len(n_pred_noise))))
linpreds_surv <- -0.1 +
    dat_sim_surv[, seq_len(n_pred_truth), drop = FALSE] %*%
    rep_len(c(0.3, -0.2), length.out = n_pred_truth)
epreds_surv <- exp(linpreds_surv)
dat_sim_surv <- as.data.frame(dat_sim_surv)
cens_surv <- runif(N_surv,
                   min = quantile(epreds_surv, probs = 0.4),
                   max = quantile(epreds_surv, probs = 0.9))

## ----weibull_data-------------------------------------------------------------
shape_weib <- 1.2
scales_weib <- epreds_surv / (gamma(1 + (1 / shape_weib)))
y_weib <- rweibull(N_surv, shape = shape_weib, scale = scales_weib)
is_event_weib <- y_weib < cens_surv
yobs_weib <- y_weib
yobs_weib[!is_event_weib] <- cens_surv[!is_event_weib]
dat_sim_weib <- data.frame(yobs = yobs_weib,
                           is_censored = 1 - is_event_weib,
                           dat_sim_surv)

## ----weibull_fit--------------------------------------------------------------
refm_fit_weib <- brms::brm(
  formula = yobs | cens(is_censored) ~ .,
  family = brms::weibull(),
  data = dat_sim_weib,
  prior = brms::prior(R2D2(mean_R2 = 1/3, prec_R2 = 3)),
  ### Only for the sake of speed (not recommended in general):
  chains = 2,
  ###
  silent = 2,
  refresh = 0
)

## ----weibull_prep_projpred----------------------------------------------------
refm_shape <- as.matrix(refm_fit_weib)[, "shape", drop = FALSE]

latent_ll_oscale_weib <- structure(function(
    ilpreds,
    dis = rep(NA, nrow(ilpreds)),
    y_oscale,
    wobs = rep(1, ncol(ilpreds)),
    cens,
    cl_ref,
    wdraws_ref = rep(1, length(cl_ref))
) {
  idxs_cens <- which(cens == 1)
  idxs_event <- setdiff(seq_along(cens), idxs_cens)
  wobs_mat <- matrix(wobs, nrow = nrow(ilpreds), ncol = ncol(ilpreds),
                     byrow = TRUE)
  refm_shape_agg <- cl_agg(refm_shape, cl = cl_ref, wdraws = wdraws_ref)
  ll_unw <- matrix(nrow = nrow(ilpreds), ncol = ncol(ilpreds))
  for (idx_cens in idxs_cens) {
    ll_unw[, idx_cens] <- pweibull(
      y_oscale[idx_cens],
      shape = refm_shape_agg,
      scale = ilpreds[, idx_cens] / gamma(1 + 1 / as.vector(refm_shape_agg)),
      lower.tail = FALSE,
      log.p = TRUE
    )
  }
  for (idx_event in idxs_event) {
    ll_unw[, idx_event] <- dweibull(
      y_oscale[idx_event],
      shape = refm_shape_agg,
      scale = ilpreds[, idx_event] / gamma(1 + 1 / as.vector(refm_shape_agg)),
      log = TRUE
    )
  }
  return(wobs_mat * ll_unw)
}, cens_var = ~ is_censored)

latent_ppd_oscale_weib <- function(
    ilpreds_resamp,
    dis_resamp = rep(NA, nrow(ilpreds_resamp)),
    wobs = rep(1, ncol(ilpreds_resamp)),
    cl_ref,
    wdraws_ref = rep(1, length(cl_ref)),
    idxs_prjdraws
) {
  warning("The draws from this `latent_ppd_oscale` function are uncensored.")
  refm_shape_agg <- cl_agg(refm_shape, cl = cl_ref, wdraws = wdraws_ref)
  refm_shape_agg_resamp <- refm_shape_agg[idxs_prjdraws, , drop = FALSE]
  ppd <- rweibull(
    prod(dim(ilpreds_resamp)),
    shape = refm_shape_agg_resamp,
    scale = ilpreds_resamp / gamma(1 + 1 / as.vector(refm_shape_agg_resamp))
  )
  ppd <- matrix(ppd, nrow = nrow(ilpreds_resamp), ncol = ncol(ilpreds_resamp))
  return(ppd)
}

refm_weib <- get_refmodel(
  refm_fit_weib,
  latent = TRUE,
  latent_ll_oscale = latent_ll_oscale_weib,
  latent_ppd_oscale = latent_ppd_oscale_weib
)

## ----weibull_cvvs-------------------------------------------------------------
cvvs_weib <- cv_varsel(
  refm_weib,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0
)

## ----weibull_plot_cvvs, fig.alt="Plot of the MLPD difference between each submodel and the reference model, with uncertainty bars, for the Weibull model with right-censored observations. The top x-axis shows the submodel size and the bottom x-axis the predictor added at each size. The first 10 predictors are the 10 truly relevant ones. The difference increases with the size and then levels off close to zero (dashed line)."----
plot(cvvs_weib, stats = "mlpd", deltas = TRUE)

## ----weibull_size-------------------------------------------------------------
size_decided_weib <- 10

## ----weibull_size_sgg---------------------------------------------------------
suggest_size(cvvs_weib, stat = "mlpd")

## ----weibull_pppc, fig.alt="Kaplan-Meier overlay plot with the time on the x-axis and the survival probability on the y-axis: the Kaplan-Meier curve of the observed times as a dark line, and the Kaplan-Meier curves of predictive draws from the Weibull submodel as light lines. The light lines are close to the dark line."----
predictors_final_weib <- head(ranking(cvvs_weib)[["fulldata"]],
                              size_decided_weib)
prj_weib <- project(refm_weib, predictor_terms = predictors_final_weib)
prj_predict_weib <- proj_predict(prj_weib)
bayesplot::ppc_km_overlay(y = dat_sim_weib$yobs, yrep = prj_predict_weib,
                          status_y = 1 - dat_sim_weib$is_censored)

## ----lognormal_data-----------------------------------------------------------
sdlog_lnorm <- 0.3
y_lnorm <- rlnorm(N_surv, meanlog = linpreds_surv, sdlog = sdlog_lnorm)
is_event_lnorm <- y_lnorm < cens_surv
yobs_lnorm <- y_lnorm
yobs_lnorm[!is_event_lnorm] <- cens_surv[!is_event_lnorm]
dat_sim_lnorm <- data.frame(yobs = yobs_lnorm,
                            is_censored = 1 - is_event_lnorm,
                            dat_sim_surv)

## ----lognormal_fit------------------------------------------------------------
refm_fit_lnorm <- brms::brm(
  formula = yobs | cens(is_censored) ~ .,
  family = brms::lognormal(),
  data = dat_sim_lnorm,
  prior = brms::prior(R2D2(mean_R2 = 1/3, prec_R2 = 3)),
  ### Only for the sake of speed (not recommended in general):
  chains = 2,
  ###
  silent = 2,
  refresh = 0
)

## ----lognormal_prep_projpred--------------------------------------------------
latent_ll_oscale_lnorm <- structure(function(
    ilpreds,
    dis = rep(NA, nrow(ilpreds)),
    y_oscale,
    wobs = rep(1, ncol(ilpreds)),
    cens,
    cl_ref,
    wdraws_ref = rep(1, length(cl_ref))
) {
  idxs_cens <- which(cens == 1)
  idxs_event <- setdiff(seq_along(cens), idxs_cens)
  wobs_mat <- matrix(wobs, nrow = nrow(ilpreds), ncol = ncol(ilpreds),
                     byrow = TRUE)
  ll_unw <- matrix(nrow = nrow(ilpreds), ncol = ncol(ilpreds))
  for (idx_cens in idxs_cens) {
    ll_unw[, idx_cens] <- plnorm(
      y_oscale[idx_cens],
      meanlog = ilpreds[, idx_cens],
      sdlog = dis,
      lower.tail = FALSE,
      log.p = TRUE
    )
  }
  for (idx_event in idxs_event) {
    ll_unw[, idx_event] <- dlnorm(
      y_oscale[idx_event],
      meanlog = ilpreds[, idx_event],
      sdlog = dis,
      log = TRUE
    )
  }
  return(wobs_mat * ll_unw)
}, cens_var = ~ is_censored)

latent_ppd_oscale_lnorm <- function(
    ilpreds_resamp,
    dis_resamp = rep(NA, nrow(ilpreds_resamp)),
    wobs = rep(1, ncol(ilpreds_resamp)),
    cl_ref,
    wdraws_ref = rep(1, length(cl_ref)),
    idxs_prjdraws
) {
  warning("The draws from this `latent_ppd_oscale` function are uncensored.")
  ppd <- rlnorm(
    prod(dim(ilpreds_resamp)),
    meanlog = ilpreds_resamp,
    sdlog = dis_resamp
  )
  ppd <- matrix(ppd, nrow = nrow(ilpreds_resamp), ncol = ncol(ilpreds_resamp))
  return(ppd)
}

refm_lnorm <- get_refmodel(
  refm_fit_lnorm,
  latent = TRUE,
  latent_ilink = function(lpreds, cl_ref,
                          wdraws_ref = rep(1, length(cl_ref))) lpreds,
  latent_ll_oscale = latent_ll_oscale_lnorm,
  latent_ppd_oscale = latent_ppd_oscale_lnorm,
  dis = as.matrix(refm_fit_lnorm)[, "sigma", drop = FALSE]
)

## ----lognormal_cvvs-----------------------------------------------------------
cvvs_lnorm <- cv_varsel(
  refm_lnorm,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0
)

## ----lognormal_plot_cvvs, fig.alt="The same kind of plot for the log-normal model with right-censored observations. The first 10 predictors are the 10 truly relevant ones. The difference increases with the size and then levels off close to zero (dashed line)."----
plot(cvvs_lnorm, stats = "mlpd", deltas = TRUE)

## ----lognormal_size-----------------------------------------------------------
size_decided_lnorm <- 10

## ----lognormal_size_sgg-------------------------------------------------------
suggest_size(cvvs_lnorm, stat = "mlpd")

## ----lognormal_pppc, fig.alt="The same Kaplan-Meier overlay plot for the log-normal submodel. The light lines are close to the dark line."----
predictors_final_lnorm <- head(ranking(cvvs_lnorm)[["fulldata"]],
                               size_decided_lnorm)
prj_lnorm <- project(refm_lnorm, predictor_terms = predictors_final_lnorm)
prj_predict_lnorm <- proj_predict(prj_lnorm)
bayesplot::ppc_km_overlay(y = dat_sim_lnorm$yobs, yrep = prj_predict_lnorm,
                          status_y = 1 - dat_sim_lnorm$is_censored)

## ----kfold_custom_nebin-------------------------------------------------------
get_refmodel_nebin <- function(fit, ...) {
  refm_prec <- as.matrix(fit)[, "reciprocal_dispersion", drop = FALSE]
  latent_ll_oscale_nebin <- function(ilpreds,
                                     dis = rep(NA, nrow(ilpreds)),
                                     y_oscale,
                                     wobs = rep(1, ncol(ilpreds)),
                                     cens,
                                     cl_ref,
                                     wdraws_ref = rep(1, length(cl_ref))) {
    y_oscale_mat <- matrix(y_oscale, nrow = nrow(ilpreds),
                           ncol = ncol(ilpreds), byrow = TRUE)
    wobs_mat <- matrix(wobs, nrow = nrow(ilpreds), ncol = ncol(ilpreds),
                       byrow = TRUE)
    refm_prec_agg <- cl_agg(refm_prec, cl = cl_ref, wdraws = wdraws_ref)
    ll_unw <- dnbinom(y_oscale_mat, size = refm_prec_agg, mu = ilpreds,
                      log = TRUE)
    return(wobs_mat * ll_unw)
  }
  latent_ppd_oscale_nebin <- function(
      ilpreds_resamp,
      dis_resamp = rep(NA, nrow(ilpreds_resamp)),
      wobs = rep(1, ncol(ilpreds_resamp)),
      cl_ref,
      wdraws_ref = rep(1, length(cl_ref)),
      idxs_prjdraws
  ) {
    refm_prec_agg <- cl_agg(refm_prec, cl = cl_ref, wdraws = wdraws_ref)
    refm_prec_agg_resamp <- refm_prec_agg[idxs_prjdraws, , drop = FALSE]
    ppd <- rnbinom(prod(dim(ilpreds_resamp)), size = refm_prec_agg_resamp,
                   mu = ilpreds_resamp)
    ppd <- matrix(ppd, nrow = nrow(ilpreds_resamp), ncol = ncol(ilpreds_resamp))
    return(ppd)
  }
  get_refmodel(fit, latent = TRUE,
               latent_ll_oscale = latent_ll_oscale_nebin,
               latent_ppd_oscale = latent_ppd_oscale_nebin, ...)
}
refm_nebin_kfold <- get_refmodel_nebin(
  refm_fit_nebin,
  cvrefbuilder = function(cvfit) {
    get_refmodel_nebin(cvfit, called_from_cvrefbuilder = TRUE)
  }
)

