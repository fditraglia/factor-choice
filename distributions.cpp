// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
colvec draw_normal(colvec mu, mat Sigma_inv){
/*-------------------------------------------------------
# RETURNS: 
#  A draw from a MV Normal(mu, Sigma_inv) distribution 
#--------------------------------------------------------
# ARGUMENTS:
#  mu           mean vector
#  Sigma_inv    precision matrix (inverse of cov matrix)
#--------------------------------------------------------
# NOTE: parameterized using precision matrix!
#-------------------------------------------------------*/
  RNGScope scope;
  int p = Sigma_inv.n_cols;
  colvec x = rnorm(p);
  mat R = chol(Sigma_inv);
  return mu + solve(trimatu(R), x);
}

// [[Rcpp::export]]
double density_normal(colvec x, colvec mu, mat Sigma_inv){
/*-------------------------------------------------------
# RETURNS: 
#  MV Normal(mu, Sigma_inv) probability density function
#--------------------------------------------------------
# ARGUMENTS:
#  x            point at which density is evaluated
#  mu           mean vector
#  Sigma_inv    precision matrix (inverse of cov matrix)
#--------------------------------------------------------
# NOTES: 
#  (1) Parameterized using precision matrix
#  (2) Intermediate steps calculated in logs for stability
#-------------------------------------------------------*/
 int p = Sigma_inv.n_cols;
 double first = -0.5 * p * log(2.0 * datum::pi);
 double val;
 double sign;
 log_det(val, sign, Sigma_inv);
 double second = 0.5 * val;
 double third = -0.5 * as_scalar(trans(x - mu) * Sigma_inv * (x - mu));
 return exp(first + second + third);
}

// [[Rcpp::export]]
mat draw_wishart(int v, mat S){
/*-------------------------------------------------------
# RETURNS: 
#  A draw from the Wishart(v, S) distribution. 
#--------------------------------------------------------
# ARGUMENTS:
#  v     degrees of freedom
#  S     scale matrix 
#--------------------------------------------------------
# NOTES: 
#  (1) Employs Bartlett Decomp. (Smith & Hocking, 1972)
#  (2) Output is identical to rwish from MCMCpack R
#      package provided the same seed is used.
#-------------------------------------------------------*/
  RNGScope scope;
  int p = S.n_rows;
  mat L = chol(S, "lower");
  mat A(p,p, fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df)); 
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  mat LA = trimatl(trimatl(L) * trimatl(A));
  return LA * LA.t();
}


// [[Rcpp::export]]
double log_mv_gamma(int p, double a){
/*-------------------------------------------------------
# RETURNS: 
#  Natural logarithm of p-dimensional MV Gamma function
#--------------------------------------------------------
# ARGUMENTS:
#  p     dimension of MV Gamma function
#  a     argument of MV Gamma function 
#--------------------------------------------------------
# NOTES: The multivariate Gamma function appears in the 
#        normalizing constant for the Wishart distribution
#-------------------------------------------------------*/
  double lgamma_sum = 0.0;
  for(int j = 1; j <= p; j++){
    lgamma_sum += R::lgammafn(a + 0.5 * (1 - j));
  }
  return 0.25 * p * (p - 1) * log(datum::pi) + lgamma_sum; 
}


// [[Rcpp::export]]
double density_wishart(mat X, int v, mat S){
/*-------------------------------------------------------
# RETURNS: 
#  Wishart(v, S) density evaluated at X
#--------------------------------------------------------
# ARGUMENTS:
#  v     degrees of freedom
#  S     scale matrix 
#--------------------------------------------------------
# NOTES: 
#-------------------------------------------------------*/
  int p = S.n_rows;
  double X_val, X_sign;
  log_det(X_val, X_sign, X);
  double term1 = 0.5 * (v - p - 1) * X_val;
  double term2 = -0.5 * trace(solve(trimatl(S), trimatl(X)));
  double term3 = -0.5 * v * p * log(2.0);
  double S_val, S_sign;
  log_det(S_val, S_sign, S);
  double term4 = -0.5 * v * S_val;
  double term5 = -1.0 * log_mv_gamma(p, 0.5 * v);
  return exp(term1 + term2 + term3 + term4 + term5);
}
