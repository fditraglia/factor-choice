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

double density_normal(colved x, colvec mu, mat Sigma_inv){
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
 double first = -p * log(2.0 * datum::pi);
 double val;
 double sign;
 double second = log_det(Sigma_inv, val, sign);
 double third = -as_scalar(trans(x - mu) * Sigma_inv * (x - mu));
 return exp((first + second + third) / 2.0);
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
  mat D = chol(S, "lower");
  mat A(p,p, fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df)); 
  }
  for(int col = 1; col < p; col++){
    for(int row = 0; row < col; row++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  return D.t() * A.t() * A * D;
}


double density_wishart(colvec v, mat S){
/*-------------------------------------------------------
# RETURNS: 
#  Wishart(v, S) probability density function
#--------------------------------------------------------
# ARGUMENTS:
#  v     degrees of freedom
#  S     scale matrix 
#--------------------------------------------------------
# NOTES: 
#-------------------------------------------------------*/
}