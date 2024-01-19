#' @useDynLib VBRobM
#' @import lme4 stats graphics
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm rmvnorm
#' @importFrom igraph graph_from_data_frame all_simple_paths
#' @title
#' Variational Bayesian Robust Mediation Analysis
#'
#' @description
#' Overall variational Bayesian robust mediation analysis for a direct acyclic graph.
#'
#' @param ... one or more `rREML` objects
#' @param n bootstrap replication times
#' @param q desired quantile for confidence interval
#'
#'
#' @returns
#' A list with components:
#' * `coefs` a data frame contains estimation, standard error, z statistic and p-value of Wald-type test for all paths.
#' * `mediations` a summary for direct, indirect and total effects of the fitted models.
#' @export
#'
#' @examples
#' res <- vbrobm(vbrob(smd_enemy~1+(1|Code), data = trophic, vi = trophic$vi_enemy, maxit = 500),
#'               vbrob(smd_herbivore~1+smd_enemy+(1|Code), data = trophic, vi = trophic$vi_herbivore, maxit = 500),
#'               vbrob(smd_plant~1+smd_herbivore+(1|Code), data = trophic, vi = trophic$vi_plant, maxit = 500))
#' res
#' plot(res)
vbrobm <- function(..., n = 1000, q = c(0.025, 0.975)){
  x <- list(...)
  res_coefs <- do.call(rbind, lapply(x, coefs))
  res_varcovs <- as.matrix(bdiag(lapply(x, covs)))
  boo <- rmvnorm(n, res_coefs$Coefficient, res_varcovs)
  preds <- res_coefs$Predictor
  resps <- res_coefs$Response
  exposure <- setdiff(preds, resps)
  outcome <- setdiff(resps, preds)
  relations <- res_coefs[, c(1, 2)]
  g <- graph_from_data_frame(relations)
  mediation_results <- NULL
  for (i in exposure) {
    for (j in outcome) {
      paths <- all_simple_paths(g, from = i, to = j)
      mediations <- NULL
      total_betas <- rep(0, n)
      for (k in paths) {
        path_betas <- rep(1, n)
        for (l in 1:(length(k)-1)) {
          path_pair <- attr(k[c(l, l+1)], "names")
          idx <- which(rowSums('dim<-'(as.matrix(relations)%in%path_pair, dim(relations))) == 2)
          path_betas <- path_betas*boo[, idx]
        }
        total_betas <- total_betas + path_betas
        if (length(k) == 2) {
          mediations <- rbind(mediations, c(i,
                                            j,
                                            NA,
                                            "Direct effect",
                                            mean(path_betas),
                                            sd(path_betas),
                                            quantile(path_betas, q)))
        } else {
          mediations <- rbind(mediations, c(i,
                                            j,
                                            paste(attr(k, "names")[2:(length(k)-1)], collapse = '&'),
                                            "Indirect effect",
                                            mean(path_betas),
                                            sd(path_betas),
                                            quantile(path_betas, q)))
        }
      }
      mediations <- rbind(mediations, c(i,
                                        j,
                                        NA,
                                        "Total effect",
                                        mean(total_betas),
                                        sd(total_betas),
                                        quantile(total_betas, q)))
      mediations <- cbind(data.frame(mediations[, 1:4]), data.frame(apply(mediations[, 5:8], 2, as.numeric)))
      colnames(mediations) <- c("Exposure", "Outcome", "Mediators", "Effect type", "Mean effect", "Sd", "Lower bound", "Upper bound")
      mediation_results <- c(mediation_results, list(mediations))
    }
  }
  result <- list(coefs = res_coefs, mediations = mediation_results)
  class(result) <- "MedRes"
  result
}

#' @title
#' Print MedRes Object
#'
#' @description
#' Generic function to print a `MedRes` object.
#'
#' @param x a `MedRes` object to print.
#'
#'
#' @export
#'
print.MedRes <- function(x) {
  cat("Path coefficients:\n")
  print(x$coefs)
  cat("\n")
  cat("Mediation analysis:\n")
  print(x$mediations)
}


#' @title
#' Show MedRes Object
#'
#' @description
#' Generic function to show a `MedRes` object.
#'
#' @param x a `MedRes` object to show.
#' @param ... further arguments passed to or from other methods.
#'
#'
#' @export
#'
showDefault.MedRes <- function(x, ...) {
  print(x = x, ...)
}

#' @title
#' DAG Plot
#' @import ggraph
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggplot2 aes theme_void arrow unit
#' @description
#' Generic function to generate the DAG plot for a `MedRes` object.
#'
#' @param x a `MedRes` object to show.
#' @param ... further arguments passed to or from other methods.
#'
#'
#' @export
#'
plot.MedRes <- function(
    x,
    ...,
    main = "Diagnosis Plot",
    xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles") {
  res <- x$coefs[, c(1:3, 6)]
  res$significance <-ifelse(res$p.value < 0.05, "*", "")
  g <- graph_from_data_frame(res, directed = TRUE)
  ggraph(g, layout = "linear", circular = T) +
    geom_edge_link(aes(label = paste(round(Coefficient, 3), significance, sep = "")),
                   arrow = arrow(length = unit(4, 'mm'))) +
    geom_node_text(aes(label = name), vjust = -0.5, hjust = 0.5, size = 5) +
    theme_void()
}
