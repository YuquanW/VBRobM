#' @useDynLib VBRobM
#' @import lme4 stats graphics RcppArmadillo
#' @importFrom Rcpp evalCpp
#' @title
#' Robust REML Proposal II Fit
#'
#' @description
#' Robust restricted maximum likelihood (REML) proposal II estimation (only support random intercept).
#'
#' @param formula an object of class "[formula][stats::formula]" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param vi an optional vector of sampling variances to be used in the fitting process.
#' @param params an optional vector of initial guess for the parameters in the model. By default the results of "[lmer][lme4::lmer]" will be used.
#' @param k1 a positive numeric indicates the tuning parameter of bisquare loss for location estimation.
#' By default 4.685 will be used corresponding to 95% asymptotic efficiency in the simple location problem.
#' @param k2 a positive numeric indicates the tuning parameter of bisquare loss for scale estimation.
#' By default 5.12 will be used corresponding to 90% asymptotic efficiency in the simple scale problem.
#' @param maxit maximum number of iterations. Default is 500.
#' @param tol a numeric indicates convergence tolerance.
#'
#'
#' @returns
#' A list with components:
#' * `coefTable` a data frame contains estimation, standard error, z statistic and p-value of Wald-type test.
#' * `varcov` the variance-covariance matrix for the coefficients.
#' * `s1` the correction factor for sampling variance.
#' * `s2` the variance of the random intercept in 2-level meta-analysis.
#' * `Std.resid` standardized residuals of the fitted model.
#' * `convergence` 0 indicates the algorithm converged; 1 indicates the algorithm failed to converge.
#' @export
#'
#' @examples
#' fit <- rREML(smd_enemy ~ 1 + (1|Code), trophic, trophic$vi_enemy)
#' fit$coefTable
rREML <- function(formula, data, vi, params, k1 = 4.685, k2 = 5.12, maxit = 500, tol = .Machine$double.eps^0.5){
  # k2 = 2.068 for 90%; k2 = 2.376 for 95%
  # k2 = 5.12 for 90% using Tukey loss; k2 = 6.1 for 95% using Tukey

  ## warm up
  n <- nrow(data)
  if (missing(vi)) {
    vi <- rep(1, n)
  }
  environment(formula) <- environment()
  if (missing(params)) {
    results0 <- lmer(formula, data, weights = 1/vi)
    p <- results0@devcomp$dims["p"]
    ## extract data

    Xi <- getME(results0, "X")
    Ji <- tcrossprod(as.matrix(getME(results0, "Z")))
    yi <- getME(results0, "y")
    beta.new <- getME(results0, "fixef")
    s.new <- c(attr(VarCorr(results0)[[1]], "stddev")^2, sigma(results0)^2)
    r0 <- residuals(results0, type = "pearson")
  } else {
    design_mat <- lFormula(eval(formula), data)
    Xi <- design_mat$X
    Ji <- crossprod(as.matrix(design_mat$reTrms$Zt))
    yi <- design_mat$fr[, 1]
    p <- ncol(Xi)
    beta.new <- params[1:p]
    s.new <- params[c(p+1, p+2)]
  }
  D <- diag(vi)

  K2 <- integrate(function(x) psi_bisquare(x, k2)^2*dnorm(x), -Inf, Inf)$value
  res <- rcpp_coord_descend(beta.new, s.new,
                            yi, Xi, Ji, vi,
                            maxit, k1, k2, K2, tol)
  beta.new <- drop(res$beta.new)
  s.new <- drop(res$s.new)
  convergence <- res$convergence
  r <- drop(res$r)
  U <- res$U
  V = Ji*s.new[1] + diag(vi*s.new[2])
  V_inv <- solve(V)
  mu_2 <- var(psi_bisquare(r, k1))
  eta_1 <- mean(dpsi_bisquare(r, k1))
  V2 <- mu_2*t(Xi)%*%V_inv%*%Xi
  V1 <- eta_1*t(Xi)%*%V_inv%*%Xi
  varcov <- solve(V1)%*%V2%*%t(solve(V1))
  beta.se <- sqrt(diag(varcov))

  variables <- colnames(Xi)
  Predictor <- if(length(variables[-1]) == 0) "Treatment" else c("Treatment", variables[-1])
  Response <- rep(all.vars(formula)[1], length(Predictor))

  result <- list(call = match.call(),
                 formula = formula,
                 coefTable = data.frame(Predictor = Predictor,
                                        Response = Response,
                                        Coefficient = beta.new,
                                        Std.Error = beta.se,
                                        z = beta.new/beta.se,
                                        p.value = 2*pnorm(abs(beta.new/beta.se), lower.tail = F), row.names = NULL),
                 varcov = varcov,
                 s1 = sqrt(s.new[1]),
                 s2 = sqrt(s.new[2]),
                 Std.resid = r,
                 convergence = convergence)

  class(result) <- "rREML"
  result
}

#' @useDynLib VBRobM
#' @import lme4 stats graphics
#' @importFrom Rcpp evalCpp
#' @title
#' Bias-corrected Robust REML Fit
#'
#' @description
#' Bias-corrected robust REML estimation for two-level linear mixed models (only support random intercept).
#'
#' @param formula an object of class "[formula][stats::formula]" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param vi an optional vector of sampling variances to be used in the fitting process.
#' @param params an optional vector of initial guess for the parameters in the model. By default the results of "[lmer][lme4::lmer]" will be used.
#' @param k1 a positive numeric indicates the tuning parameter of bisquare loss for location estimation.
#' Defaults to 4.685 corresponding to 95% asymptotic efficiency in the simple location problem.
#' @param k2 a positive numeric indicates the tuning parameter of bisquare loss for scale estimation.
#' Defaults to 5.12 corresponding to 90% asymptotic efficiency in the simple scale problem.
#' @param maxit maximum number of iterations. Defaults to 500.
#' @param tol a numeric indicates convergence tolerance.
#'
#'
#' @returns
#' A list with components:
#' * `coefTable` a data frame contains estimation, standard error, z statistic and p-value of Wald-type test.
#' * `varcov` the variance-covariance matrix for the coefficients.
#' * `s1` the correction factor for sampling variance.
#' * `s2` the variance of the random intercept in 2-level meta-analysis.
#' * `Std.resid` standardized residuals of the fitted model.
#' * `convergence` 0 indicates the algorithm converged; 1 indicates the algorithm failed to converge.
#' @export
#'
#' @examples
#' fit <- bcrREML(smd_enemy ~ 1 + (1|Code), trophic, trophic$vi_enemy)
#' fit$coefTable
bcrREML <- function(formula, data, vi, params, k1 = 4.685, k2 = 5.12, maxit = 500, tol = .Machine$double.eps^0.5){
  # k2 = 2.068 for 90%; k2 = 2.376 for 95%
  # k2 = 5.12 for 90% using Tukey loss; k2 = 6.1 for 95% using Tukey

  ## warm up
  n <- nrow(data)
  if (missing(vi)) {
    vi <- rep(1, n)
  }
  environment(formula) <- environment()
  if (missing(params)) {
    results0 <- lmer(formula, data, weights = 1/vi)
    p <- results0@devcomp$dims["p"]
    ## extract data

    Xi <- getME(results0, "X")
    Ji <- tcrossprod(as.matrix(getME(results0, "Z")))
    yi <- getME(results0, "y")
    beta.new <- getME(results0, "fixef")
    s.new <- c(attr(VarCorr(results0)[[1]], "stddev")^2, sigma(results0)^2)
    r0 <- residuals(results0, type = "pearson")
  } else {
    design_mat <- lFormula(eval(formula), data)
    Xi <- design_mat$X
    Ji <- crossprod(as.matrix(design_mat$reTrms$Zt))
    yi <- design_mat$fr[, 1]
    p <- ncol(Xi)
    beta.new <- params[1:p]
    s.new <- params[c(p+1, p+2)]
  }
  D <- diag(vi)

  K2 <- integrate(function(x) psi_bisquare(x, k2)^2*dnorm(x), -Inf, Inf)$value
  res <- rcpp_coord_descend(beta.new, s.new,
                            yi, Xi, Ji, vi,
                            maxit, k1, k2, K2, tol)
  beta.new <- drop(res$beta.new)
  s.new <- drop(res$s.new)
  convergence <- res$convergence
  U <- res$U
  if (missing(params)) {
    beta0 <- results0@beta
  } else {
    beta0 <- lmer(formula, data, weights = 1/vi)@beta
  }
  z <- U%*%(yi-Xi%*%beta0)
  V <- Ji*s.new[1] + D*s.new[2]
  P <- solve(V) - solve(V, Xi)%*%solve(t(Xi)%*%solve(V, Xi))%*%t(solve(V, Xi))
  V_inv <- solve(V)
  dP1 <- -P%*%Ji%*%P
  dP2 <- -P%*%D%*%P
  score_beta <- drop(t(Xi)%*%U%*%psi_bisquare(z, k1))
  score_s1 <- (drop(t(psi_bisquare(z, k2))%*%U%*%Ji%*%U%*%psi_bisquare(z, k2))-
                 sum(diag(K2*P%*%Ji)))
  score_s2 <- (drop(t(psi_bisquare(z, k2))%*%U%*%D%*%U%*%psi_bisquare(z, k2))-
                 sum(diag(K2*P%*%D)))
  score <- c(score_beta, score_s1, score_s2)
  grad_betabeta <- -t(Xi)%*%U%*%diag(drop(dpsi_bisquare(z, k1)))%*%U%*%Xi
  grad_betas1 <- -0.5*drop(t(Xi)%*%U%*%Ji%*%V_inv%*%psi_bisquare(z, k1))-
    0.5*drop(t(Xi)%*%U%*%diag(drop(dpsi_bisquare(z, k1)))%*%U%*%Ji%*%U%*%z)
  grad_betas2 <- -0.5*drop(t(Xi)%*%U%*%D%*%V_inv%*%psi_bisquare(z, k1))-
    0.5*drop(t(Xi)%*%U%*%diag(drop(dpsi_bisquare(z, k1)))%*%U%*%D%*%U%*%z)
  grad_s1beta <- -2*drop(t(Xi)%*%U%*%diag(drop(dpsi_bisquare(z, k2)))%*%U%*%Ji%*%U%*%psi_bisquare(z, k2))
  grad_s1s1 <- -drop(t(z)%*%U%*%Ji%*%U%*%diag(drop(dpsi_bisquare(z, k2)))%*%U%*%Ji%*%U%*%psi_bisquare(z, k2))-
    drop(t(psi_bisquare(z, k2))%*%U%*%Ji%*%V_inv%*%Ji%*%U%*%psi_bisquare(z, k2))-
    sum(diag(K2*dP1%*%Ji))
  grad_s1s2 <- -drop(t(z)%*%U%*%D%*%U%*%diag(drop(dpsi_bisquare(z, k2)))%*%U%*%Ji%*%U%*%psi_bisquare(z, k2))-
    drop(t(psi_bisquare(z, k2))%*%U%*%D%*%V_inv%*%Ji%*%U%*%psi_bisquare(z, k2))-
    sum(diag(K2*dP2%*%Ji))
  grad_s2beta <- -2*drop(t(Xi)%*%U%*%diag(drop(dpsi_bisquare(z, k2)))%*%U%*%D%*%U%*%psi_bisquare(z, k2))
  grad_s2s1 <- -drop(t(z)%*%U%*%Ji%*%U%*%diag(drop(dpsi_bisquare(z, k2)))%*%U%*%D%*%U%*%psi_bisquare(z, k2))-
    drop(t(psi_bisquare(z, k2))%*%U%*%Ji%*%V_inv%*%D%*%U%*%psi_bisquare(z, k2))-
    sum(diag(K2*dP1%*%D))
  grad_s2s2 <- -drop(t(z)%*%U%*%D%*%U%*%diag(drop(dpsi_bisquare(z, k2)))%*%U%*%D%*%U%*%psi_bisquare(z, k2))-
    drop(t(psi_bisquare(z, k2))%*%U%*%D%*%V_inv%*%D%*%U%*%psi_bisquare(z, k2))-
    sum(diag(K2*dP2%*%D))
  grad <- cbind(rbind(grad_betabeta, grad_s1beta, grad_s2beta),
                c(grad_betas1, grad_s1s1, grad_s2s1),
                c(grad_betas2, grad_s1s2, grad_s2s2))
  dlt <- solve(grad, score)
  beta.new <- beta.new + dlt[1:p]
  s.new <- s.new + dlt[c(p+1, p+2)]
  s.new <- pmax(s.new, 0)
  V <- Ji*s.new[1] + D*s.new[2]
  V.eigen <- eigen(V)
  U <- V.eigen$vectors%*%diag(1/sqrt(V.eigen$values))%*%t(V.eigen$vectors)
  r <- drop(U%*%(yi-Xi%*%beta.new))
  mu_2 <- var(psi_bisquare(r, k1))
  eta_1 <- mean(dpsi_bisquare(r, k1))
  V2 <- mu_2*t(Xi)%*%V_inv%*%Xi
  V1 <- eta_1*t(Xi)%*%V_inv%*%Xi
  varcov <- solve(V1)%*%V2%*%t(solve(V1))
  beta.se <- sqrt(diag(varcov))

  variables <- colnames(Xi)
  Predictor <- if(length(variables[-1]) == 0) "Treatment" else c("Treatment", variables[-1])
  Response <- rep(all.vars(formula)[1], length(Predictor))

  result <- list(call = match.call(),
                 formula = formula,
                 coefTable = data.frame(Predictor = Predictor,
                                        Response = Response,
                                        Coefficient = beta.new,
                                        Std.Error = beta.se,
                                        z = beta.new/beta.se,
                                        p.value = 2*pnorm(abs(beta.new/beta.se), lower.tail = F), row.names = NULL),
                 varcov = varcov,
                 s1 = sqrt(s.new[1]),
                 s2 = sqrt(s.new[2]),
                 Std.resid = r,
                 convergence = convergence)

  class(result) <- "rREML"
  result
}

#' @useDynLib VBRobM
#' @import lme4 stats graphics
#' @importFrom Rcpp evalCpp
#' @title
#' Variational Bayesian Robust Estimating Equation Fit
#'
#' @description
#' Variational Bayesian robust estimation for two-level linear mixed models (only support random intercept).
#'
#' @param formula an object of class "[formula][stats::formula]" (or one that can be coerced to that class):
#' a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param vi an optional vector of sampling variances to be used in the fitting process.
#' @param params an optional vector of initial guess for the parameters in the model. By default the results of "[lmer][lme4::lmer]" will be used.
#' @param k1 a positive numeric indicates the tuning parameter of bisquare loss for location estimation.
#' Defaults to 4.685 corresponding to 95% asymptotic efficiency in the simple location problem.
#' @param k2 a positive numeric indicates the tuning parameter of bisquare loss for scale estimation.
#' Defaults to 5.12 corresponding to 90% asymptotic efficiency in the simple scale problem.
#' @param maxit maximum number of iterations. Defaults to 500.
#' @param tol a numeric indicates convergence tolerance.
#'
#'
#' @returns
#' A list with components:
#' * `coefTable` a data frame contains estimation, standard error, z statistic and p-value of Wald-type test.
#' * `varcov` the variance-covariance matrix for the coefficients.
#' * `s1` the correction factor for sampling variance.
#' * `s2` the variance of the random intercept in 2-level meta-analysis.
#' * `Std.resid` standardized residuals of the fitted model.
#' * `convergence` 0 indicates the algorithm converged; 1 indicates the algorithm failed to converge.
#' @export
#'
#' @examples
#' fit <- vbrob(smd_enemy ~ 1 + (1|Code), trophic, trophic$vi_enemy)
#' fit$coefTable
vbrob <- function(formula, data, vi, params, k1 = 4.685, k2 = 5.12, maxit = 500, tol = .Machine$double.eps^0.5){
  # k2 = 2.068 for 90%; k2 = 2.376 for 95%
  # k2 = 5.12 for 90% using Tukey loss; k2 = 6.1 for 95% using Tukey

  n <- nrow(data)
  if (missing(vi)) {
    vi <- rep(1, n)
  }
  environment(formula) <- environment()
  if (missing(params)) {
    results0 <- lmer(formula, data, weights = 1/vi)
    p <- results0@devcomp$dims["p"]
    ## extract data

    Xi <- getME(results0, "X")
    Ji <- tcrossprod(as.matrix(getME(results0, "Z")))
    yi <- getME(results0, "y")
    beta.new <- getME(results0, "fixef")
    s.new <- c(attr(VarCorr(results0)[[1]], "stddev")^2, sigma(results0)^2)
    r0 <- residuals(results0, type = "pearson")
  } else {
    design_mat <- lFormula(eval(formula), data)
    Xi <- design_mat$X
    Ji <- crossprod(as.matrix(design_mat$reTrms$Zt))
    yi <- design_mat$fr[, 1]
    p <- ncol(Xi)
    beta.new <- params[1:p]
    s.new <- params[c(p+1, p+2)]
  }
  D <- diag(vi)

  K2 <- integrate(function(x) psi_bisquare(x, k2)^2*dnorm(x), -Inf, Inf)$value
  res <- rcpp_coord_descend(beta.new, s.new,
                            yi, Xi, Ji, vi,
                            maxit, k1, k2, K2, tol)
  beta.new <- drop(res$beta.new)
  s.new <- drop(res$s.new)
  U <- res$U
  # if (missing(params)) {
  #   beta0 <- results0@beta
  # } else {
  #   beta0 <- lmer(formula, data, weights = 1/vi)@beta
  # }
  dpmm <- nig_dpmm(drop(U%*%(yi-Xi%*%beta.new)), tau = 100, s = 1, t = 1, alpha = 1, K = 20, maxiter = maxit)$mixpdf
  K1 <- integrate(function(x)psi_bisquare(x, k1)*mapply(dpmm, x), -Inf, Inf, rel.tol = tol)$value
  K2 <- integrate(function(x)psi_bisquare(x, k2)*mapply(dpmm, x), -Inf, Inf, rel.tol = tol)$value
  K3 <- integrate(function(x)psi_bisquare(x, k2)^2*mapply(dpmm, x), -Inf, Inf, rel.tol = tol)$value - K2^2
  res <- rcpp_coord_descend_stg2(beta.new, s.new,
                                 yi, Xi, Ji, vi,
                                 maxit, k1, k2, K1, K2, K3, tol)
  beta.new <- drop(res$beta.new)
  s.new <- drop(res$s.new)
  convergence <- res$convergence
  r <- res$r[, 1]
  U <- res$U
  V_inv <- U%*%U
  mu_2 <- var(psi_bisquare(r, k1))
  eta_1 <- mean(dpsi_bisquare(r, k1))
  V2 <- mu_2*t(Xi)%*%V_inv%*%Xi
  V1 <- eta_1*t(Xi)%*%V_inv%*%Xi
  varcov <- solve(V1)%*%V2%*%t(solve(V1))
  beta.se <- sqrt(diag(varcov))

  variables <- colnames(Xi)
  Predictor <- if(length(variables[-1]) == 0) "Treatment" else c("Treatment", variables[-1])
  Response <- rep(all.vars(formula)[1], length(Predictor))

  result <- list(call = match.call(),
                 formula = formula,
                 coefTable = data.frame(Predictor = Predictor,
                                        Response = Response,
                                        Coefficient = beta.new,
                                        Std.Error = beta.se,
                                        z = beta.new/beta.se,
                                        p.value = 2*pnorm(abs(beta.new/beta.se), lower.tail = F), row.names = NULL),
                 varcov = varcov,
                 s1 = sqrt(s.new[1]),
                 s2 = sqrt(s.new[2]),
                 Std.resid = r,
                 convergence = convergence)

  class(result) <- "rREML"
  result
}

#' @title
#' Extract Fixed Coefficients
#'
#' @description
#' A function for extracting fixed coefficients in the model.
#'
#' @param x a `rREML` object.
#'
#' @examples
#' x <- rREML(smd_enemy ~ 1 + (1 | Code), trophic, trophic$vi_enemy)
#' coefs(x)
#'
#' @export
#'
coefs <- function(x) {
  x$coefTable
}

#' @title
#' Extract Variance-Covariance Matrix
#'
#' @description
#' A function for extracting variance-covariance matrix for fixed coefficients in the model.
#'
#' @param x a `rREML` object.
#'
#' @examples
#' x <- rREML(smd_enemy ~ 1 + (1 | Code), trophic, trophic$vi_enemy)
#' covs(x)
#'
#' @export
#'
covs <- function(x) {
  x$varcov
}

#' @title
#' Model Predictions
#'
#' @description
#' A generic function for predicting within-level value of response.
#'
#' @param x a `rREML` object.
#' @param newdata a data frame contains same variables as fitted object x.
#'
#' @examples
#' x <- rREML(smd_enemy ~ 1 + (1 | Code), trophic, trophic$vi_enemy)
#' pred(x, trophic)
#'
#' @export
#'
pred <- function(x, newdata) {
  formula <- x$formula
  design_mat <- lFormula(eval(formula), newdata)
  X <- design_mat$X
  pred <- drop(X%*%x$coefTable$Coefficient)
  pred
}

#' @title
#' Print rREML Object
#'
#' @description
#' Generic function to print a `rREML` object.
#'
#' @param x a `rREML` object to print.
#'
#' @examples
#' x <- rREML(smd_enemy ~ 1 + (1 | Code), trophic, trophic$vi_enemy)
#' print(x)
#'
#' @export
#'
print.rREML <- function(x, ...) {
  cat("Call:\n")
  cat(deparse(x$call), "\n", "\n")
  cat("Coefficients:\n")
  x$coefTable[, c(-1, -2)] <- round(x$coefTable[, c(-1, -2)], digits = 3)
  print(format(x$coefTable), nsmall = 3)
  cat("\n")
  cat("Between-level Std.Error:\n")
  cat(x$s1, "\n")
  cat("\n")
  cat("Convergence\n")
  cat(x$convergence)
}


#' @title
#' Show rREML Object
#'
#' @description
#' Generic function to show a `rREML` object.
#'
#' @param x a `rREML` object to show.
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#' rREML(smd_enemy ~ 1 + (1 | Code), trophic, trophic$vi_enemy)
#'
#' @export
#'
showDefault.rREML <- function(x, ...) {
  print(x = x, ...)
}

#' @title
#' Diagnosis Plot
#'
#' @description
#' Generic function to generate the diagnosis plot for a `rREML` object.
#'
#' @param x a `rREML` object to show.
#' @param ... further arguments passed to or from other methods.
#' @param main the main title for the plot.
#' @param xlab the label for x axis.
#' @param ylab the label for y axis.
#'
#' @examples
#' x <- rREML(smd_enemy ~ 1 + (1 | Code), trophic, trophic$vi_enemy)
#' plot(x)
#'
#' @export
#'
plot.rREML <- function(
    x,
    ...,
    main = "Diagnosis Plot",
    xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles") {
  qqnorm(x$Std.resid,
         ...,
         pch = 20,
         main = main,
         xlab = xlab,
         ylab = ylab)
  abline(a = 0, b = 1, col = "red")
}


