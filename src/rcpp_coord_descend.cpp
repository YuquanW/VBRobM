#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;
using namespace arma;


// [[Rcpp::export]]
arma::vec rcpp_weight_bisquare(arma::vec x, double k) {

    //int p = x.n_elem;
    //arma::vec w = arma::min(ones(p), k/abs(x));
    arma::vec w = pow(1 - pow(x/k, 2), 2)%(abs(x) <= k);
    return w;
}

// [[Rcpp::export]]
arma::vec rcpp_psi_bisquare(arma::vec x, double k) {

    //arma::vec psi = x%(abs(x) <= k) + k*sign(x)%(abs(x) > k);
    arma::vec psi = x%pow(1 - pow(x/k, 2), 2)%(abs(x) <= k);
    return psi;
}


// [[Rcpp::export]]
List rcpp_coord_descend(arma::vec beta0, arma::vec s0,
                        arma::vec Yi, arma::mat Xi, arma::mat Ji, arma::vec vi,
                        int n, double k1, double k2, double K2, double tol) {
    mat D = diagmat(vi);
    int i;
    for (int i = 0; i < n; i++) {
      vec s0_old = s0;
      mat V = Ji*s0_old(0) + D*s0_old(1);
      mat U = sqrtmat_sympd(V).i();
      for (int j = 0; j < n; j++){
        vec beta0_old = beta0;
        vec r = U*(Yi - Xi*beta0_old);
        mat W = diagmat(rcpp_weight_bisquare(r, k1));
        beta0 = solve(Xi.t()*U*W*U*Xi, Xi.t()*U*W*U*Yi);
        if (norm(beta0 - beta0_old, 2) <= tol) {
          break;
        }
      }
      vec r = U*(Yi - Xi*beta0);
      mat vinv_x = solve(V, Xi);
      mat vinv_d = solve(V, D);
      mat vinv_j = solve(V, Ji);
      mat P = V.i() - vinv_x*(Xi.t()*vinv_x).i()*vinv_x.t();
      double L00 = trace(K2*P*Ji*vinv_j);
      double L01 = trace(K2*P*Ji*vinv_d);
      double L10 = trace(K2*P*D*vinv_j);
      double L11 = trace(K2*P*D*vinv_d);
      vec R0 = rcpp_psi_bisquare(r, k2).t()*U*Ji*U*rcpp_psi_bisquare(r, k2);
      vec R1 = rcpp_psi_bisquare(r, k2).t()*U*D*U*rcpp_psi_bisquare(r, k2);
      mat L = {{L00, L01},
               {L10, L11}};
      vec R = {R0(0), R1(0)};
      s0 = max(solve(L, R), zeros(2));
      if (norm(s0 - s0_old, 2) <= tol) {
        break;
      }
    }
    int convergence = 1*(i == n);
    mat V = Ji*s0(0) + D*s0(1);
    mat U = sqrtmat_sympd(V).i();
    vec r = U*(Yi - Xi*beta0);
    return List::create(_["beta.new"] = beta0,
                        _["s.new"] = s0,
                        _["convergence"] = convergence,
                        _["r"] = r,
                        _["U"] = U);
}

// [[Rcpp::export]]
List rcpp_coord_descend_stg2(arma::vec beta0, arma::vec s0,
                             arma::vec Yi, arma::mat Xi, arma::mat Ji, arma::vec vi,
                             int n, double k1, double k2, double K1, double K2, double K3,
                             double tol) {
  mat D = diagmat(vi);
  int i;
  for (int i = 0; i < n; i++) {
    vec s0_old = s0;
    mat V = Ji*s0_old(0) + D*s0_old(1);
    mat U = sqrtmat_sympd(V).i();
    for (int j = 0; j < n; j++){
      vec beta0_old = beta0;
      vec r = U*(Yi - Xi*beta0_old);
      mat W = diagmat(rcpp_weight_bisquare(r, k1));
      vec mu(Xi.n_rows);
      mu.fill(K1);
      beta0 = solve(Xi.t()*U*W*U*Xi, Xi.t()*U*W*U*Yi-Xi.t()*U*mu);
      if (norm(beta0 - beta0_old, 2) <= tol) {
        break;
      }
    }
    vec r = U*(Yi - Xi*beta0);
    mat vinv_x = solve(V, Xi);
    mat vinv_d = solve(V, D);
    mat vinv_j = solve(V, Ji);
    mat P = V.i() - vinv_x*(Xi.t()*vinv_x).i()*vinv_x.t();
    double L00 = trace(K3*P*Ji*vinv_j);
    double L01 = trace(K3*P*Ji*vinv_d);
    double L10 = trace(K3*P*D*vinv_j);
    double L11 = trace(K3*P*D*vinv_d);
    vec mu(Xi.n_rows);
    mu.fill(K2);
    vec R0 = rcpp_psi_bisquare(r, k2).t()*U*Ji*U*rcpp_psi_bisquare(r, k2) - mu.t()*U*Ji*U*mu;
    vec R1 = rcpp_psi_bisquare(r, k2).t()*U*D*U*rcpp_psi_bisquare(r, k2) - mu.t()*U*D*U*mu;
    mat L = {{L00, L01},
             {L10, L11}};
    vec R = {R0(0), R1(0)};
    s0 = max(solve(L, R), zeros(2));
    if (norm(s0 - s0_old, 2) <= tol) {
      break;
    }
  }
  int convergence = 1*(i == n);
  mat V = Ji*s0(0) + D*s0(1);
  mat U = sqrtmat_sympd(V).i();
  vec r = U*(Yi - Xi*beta0);
  return List::create(_["beta.new"] = beta0,
                      _["s.new"] = s0,
                      _["convergence"] = convergence,
                      _["r"] = r,
                      _["U"] = U);
}
