#' @import stats
#' @export
nig_dpmm <- function(data, K = 20, tau = 25, s = 1, t = 0.04, alpha = 1, maxiter = 100, tol = 1e-5)
{
  n <- length(data)
  p_k <- matrix(1/K, n, K)
  n_k <- colSums(p_k)
  n_k_greater <- replicate(K, n_k)
  n_k_greater[upper.tri(n_k_greater, diag = T)] <- 0
  n_k_greater <- colSums(n_k_greater)
  mu_k <- colSums(p_k*data)/(1/tau+n_k)
  tau_k <- (1/tau+n_k)^(-1)
  s_k <- s+n_k/2
  T_k <- t+0.5*(drop(t(data^2)%*%p_k)-1/tau_k*mu_k^2)
  param1 <- c(mu_k, tau_k, s_k, T_k)
  for (i in 1:maxiter) {
    param0 <- param1
    n_k_greater_less <- replicate(K, 1/(alpha+n_k_greater))
    n_k_greater_less[lower.tri(n_k_greater_less, diag = T)] <- 0
    n_k_greater_less <- colSums(n_k_greater_less)
    p_k <- data%*%t(mu_k*s_k/T_k)-0.5*data^2%*%t(s_k/T_k)
    p_k <- t(t(p_k)-0.5*(log(T_k)-digamma(s_k)+tau_k+mu_k^2*s_k/T_k)+digamma(n_k+1)-n_k_greater_less)
    p_k[, K] <- p_k[, K]-digamma(n_k[K]+1)+digamma(n_k[K]+alpha)+1/(alpha+n_k_greater[K-1])-1/(1+n_k[K])
    p_k <- p_k - apply(p_k, 1, max)
    p_k <- exp(p_k)/rowSums(exp(p_k))
    n_k <- colSums(p_k)
    n_k_greater <- replicate(K, n_k)
    n_k_greater[upper.tri(n_k_greater, diag = T)] <- 0
    n_k_greater <- colSums(n_k_greater)
    mu_k <- colSums(p_k*data)/(1/tau+n_k)
    tau_k <- (1/tau+n_k)^(-1)
    s_k <- s+n_k/2
    T_k <- t+0.5*(drop(t(data^2)%*%p_k)-1/tau_k*mu_k^2)
    param1 <- c(mu_k, tau_k, s_k, T_k)
    if (norm(param1-param0, "2") < tol) {
      break
    }
  }
  n_k_greater <- c(n, n_k_greater)
  w_k <- rep(0, K)
  for (k in 1:K) {
    if (k == 1) {
      w_k[k] <- (n_k[k]+1)/(n_k_greater[k]+alpha+1)
    } else if (k>1 & k<K) {
      w_k[k] <- (n_k[k]+1)/(n_k_greater[k]+alpha+1)*prod(1-(n_k[1:(k-1)]+1)/(n_k_greater[1:(k-1)]+alpha+1))
    } else if (k == K) {
      w_k[k] <- (n_k[k]+1)/(n_k[k]+2)*prod(1-(n_k[1:(k-1)]+1)/(n_k_greater[1:(k-1)]+alpha+1))
    }
  }
  mixpdf <- function(x)
  {
    sum(w_k*dt((x-mu_k)/(sqrt(T_k/s_k*(tau_k+1))), df = 2*s_k)/(sqrt(T_k/s_k*(tau_k+1))))
  }
  return(list(mixpdf = mixpdf, w_k = w_k, prob = p_k, N_k = n_k, mu = mu_k, tau = tau_k, s = s_k, t = T_k))
}
