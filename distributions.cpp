// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
colvec draw_normal(colvec mu, mat Sigma){
/*-------------------------------------------------------
# RETURNS: 
#  A draw from a MV Normal(mu, Sigma) distribution 
#--------------------------------------------------------
# ARGUMENTS:
#  mu       mean vector
#  Sigma    covariance matrix
#--------------------------------------------------------
# NOTES: 
#  (1) Uses eigen-decomposition (not Choleski) to match 
#      mvrnorm from the MASS package in R. 
#  (2) Output will *not* match mvrnorm even though the 
#      same standard normal draws are used: armadillo
#      and R use different conventions for eigen-decomp.
#-------------------------------------------------------*/
  RNGScope scope;
  int p = Sigma.n_cols;
  colvec x = rnorm(p);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  mat A = eigvec * diagmat(sqrt(max(eigval, zeros(p))));
  return mu + A * x;
}

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
// [[Rcpp::export]]
mat draw_wishart(int v, mat S){
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


double density_normal(colvec mu, mat Sigma_inv){
/*-------------------------------------------------------
# RETURNS: 
#  MV Normal(mu, Sigma_inv) probability density function
#--------------------------------------------------------
# ARGUMENTS:
#  mu           mean vector
#  Sigma_inv    precision matrix (inverse of cov matrix)
#--------------------------------------------------------
# NOTES: 
#-------------------------------------------------------*/
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