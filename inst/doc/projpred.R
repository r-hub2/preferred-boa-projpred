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

## ----data---------------------------------------------------------------------
data("df_gaussian", package = "projpred")
dat_gauss <- data.frame(y = df_gaussian$y, df_gaussian$x)

## ----rstanarm_attach, message=FALSE-------------------------------------------
library(rstanarm)

## ----rh-----------------------------------------------------------------------
# Number of regression coefficients:
( D <- sum(grepl("^X", names(dat_gauss))) )
# Prior guess for the number of relevant (i.e., non-zero) regression
# coefficients:
p0 <- 5
# Number of observations:
N <- nrow(dat_gauss)
# Hyperprior scale for tau, the global shrinkage parameter (note that for the
# Gaussian family, 'rstanarm' will automatically scale this by the residual
# standard deviation):
tau0 <- p0 / (D - p0) * 1 / sqrt(N)

## ----ref_fit------------------------------------------------------------------
# Set this manually if desired:
ncores <- parallel::detectCores(logical = FALSE)
### Only for technical reasons in this vignette (you can omit this when running
### the code yourself):
ncores <- min(ncores, 2L)
###
options(mc.cores = ncores)
set.seed(50780)
refm_fit <- stan_glm(
  y ~ .,
  family = gaussian(),
  data = dat_gauss,
  prior = hs(global_scale = tau0),
  ### Only for the sake of speed (not recommended in general):
  chains = 2, iter = 1000,
  ### In general, use defaults or explicitly define
  ### chains = 4, iter = 4000,
  ###
  refresh = 0
)

## ----projpred_attach, message=FALSE-------------------------------------------
library(projpred)

## ----refmodel_create----------------------------------------------------------
refm_obj <- get_refmodel(refm_fit)

## ----cv_varsel, message=FALSE-------------------------------------------------
# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
# Note: Parallel processing is disabled during package building to avoid issues
use_parallel <- FALSE  # Set to TRUE for actual parallel processing
if (use_parallel) {
  future::plan(future.mirai::mirai_multisession, workers = ncores)
}
cvvs <- cv_varsel(
  refm_obj,
  parallel = use_parallel,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0
  ###
)
# Tear down the CV parallelization setup:
if (use_parallel) {
  future::plan(future::sequential)
}

## ----plot_vsel, fig.alt="Plot of the MLPD difference between each submodel and the reference model against the submodel size from 0 to 12, with uncertainty bars. Below each size, the x-axis shows the predictor that the full-data search added and the proportion of CV folds that put this predictor at the same position. A dashed horizontal line at zero marks the reference model. The difference reaches zero at size 7 and stays slightly above zero for larger sizes."----
ggplot2::theme_set(bayesplot::theme_default(base_family = "sans"))
options(projpred.plot_vsel_size_position = "primary_x_bottom")
options(projpred.plot_vsel_text_angle = 0)
options(projpred.plot_vsel_show_cv_proportions = TRUE)
plot(cvvs, stats = "mlpd", deltas = TRUE, nterms_max = 12)

## ----size_man-----------------------------------------------------------------
size_decided <- 7

## ----size_sgg-----------------------------------------------------------------
suggest_size(cvvs, stat = "mlpd")

## ----smmry_vsel---------------------------------------------------------------
smmry <- summary(cvvs,
                 stats = "mlpd",
                 type = c("mean", "lower", "upper"),
                 deltas = TRUE,
                 nterms_max = 12)
print(smmry, digits = 2)

## ----perf_smmry---------------------------------------------------------------
perf <- performances(smmry)
str(perf)

## ----ranking------------------------------------------------------------------
rk <- ranking(cvvs, nterms_max = 12)

## ----cv_proportions-----------------------------------------------------------
( pr_rk <- cv_proportions(rk) )

## ----ranking_fulldata---------------------------------------------------------
rk[["fulldata"]]

## ----plot_cv_proportions, fig.alt="Heat map of the ranking proportions. The y-axis shows the predictors in the order of the full-data predictor ranking and the x-axis the ranking positions 1 to 12. Each tile gives the proportion of CV folds that have the predictor at that position."----
plot(pr_rk)

## ----predictors_final---------------------------------------------------------
( predictors_final <- head(rk[["fulldata"]], size_decided) )

## ----plot_cv_proportions_cumul, fig.alt="Heat map of the cumulated ranking proportions. For each predictor (y-axis) and each size j from 1 to 12 (x-axis), the tile gives the proportion of CV folds that have the predictor among their first j predictors. At sizes 7 and 8, all tiles of the first 7 or 8 predictors are at 100 %; at size 9, fewer folds include the 9th predictor."----
plot(cv_proportions(rk, cumulate = TRUE))

## ----project------------------------------------------------------------------
prj <- project(
  refm_obj,
  predictor_terms = predictors_final,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = 0
  ###
)

## ----as_matrix_prj------------------------------------------------------------
prj_mat <- as.matrix(prj)

## ----posterior_attach, message=FALSE------------------------------------------
library(posterior)

## ----smmry_prj----------------------------------------------------------------
prj_drws <- as_draws_matrix(prj_mat)
prj_smmry <- summarize_draws(
  prj_drws,
  "median", "mad", function(x) quantile(x, probs = c(0.025, 0.975))
)
# Coerce to a `data.frame` because some pkgdown versions don't print the
# tibble correctly:
prj_smmry <- as.data.frame(prj_smmry)
print(prj_smmry, digits = 1)

## ----bayesplot_attach, message=FALSE------------------------------------------
library(bayesplot)

## ----bayesplot_prj, fig.alt="Interval plot of the projected posterior, with one row for the intercept, each of the 7 selected predictors and sigma, and the parameter value on the x-axis. Each row has a point for the median, a thick line for the 50 % interval and a thin line for the 90 % interval."----
mcmc_intervals(prj_mat) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))

## ----bayesplot_ref, fig.alt="The same interval plot for the posterior of the same parameters in the reference model. The intervals differ only slightly from those of the projected posterior."----
refm_mat <- as.matrix(refm_fit)
mcmc_intervals(refm_mat, pars = colnames(prj_mat)) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))

## ----data_new-----------------------------------------------------------------
( dat_gauss_new <- setNames(
  as.data.frame(replicate(length(predictors_final), c(-1, 0, 1))),
  predictors_final
) )

## ----proj_linpred-------------------------------------------------------------
prj_linpred <- proj_linpred(prj, newdata = dat_gauss_new, integrated = TRUE)
cbind(dat_gauss_new, linpred = as.vector(prj_linpred[["pred"]]))

## ----proj_predict, fig.alt="Density overlay plot with the response value on the x-axis: the density of the observed response as a dark line, and the densities of predictive draws from the submodel as many light lines. The light lines are close to the dark line."----
prj_predict <- proj_predict(prj)
# Using the 'bayesplot' package:
ppc_dens_overlay(y = dat_gauss$y, yrep = prj_predict)

## ----ref_fit_mlvl, eval=FALSE-------------------------------------------------
# data("VerbAgg", package = "lme4")
# refm_fit <- stan_glmer(
#   r2 ~ btype + situ + mode + (btype + situ + mode | id),
#   family = binomial(),
#   data = VerbAgg,
#   QR = TRUE, refresh = 0
# )

## ----ref_fit_addv, eval=FALSE-------------------------------------------------
# data("lasrosas.corn", package = "agridat")
# # Convert `year` to a `factor` (this could also be solved by using
# # `factor(year)` in the formula, but we avoid that here to put more emphasis on
# # the demonstration of the smooth term):
# lasrosas.corn$year <- as.factor(lasrosas.corn$year)
# refm_fit <- stan_gamm4(
#   yield ~ year + topo + t2(nitro, bv),
#   family = gaussian(),
#   data = lasrosas.corn,
#   QR = TRUE, refresh = 0
# )

## ----ref_fit_addv_mlvl, eval=FALSE--------------------------------------------
# data("gumpertz.pepper", package = "agridat")
# refm_fit <- stan_gamm4(
#   disease ~ field + leaf + s(water),
#   random = ~ (1 | row) + (1 | quadrat),
#   family = binomial(),
#   data = gumpertz.pepper,
#   QR = TRUE, refresh = 0
# )

