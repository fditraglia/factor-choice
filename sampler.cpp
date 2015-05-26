// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "distributions.h"
#include "matrix_tools.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
colvec testy(colvec m){
  int p = m.n_elem;
  mat precision = eye(p, p);
  colvec out = draw_normal(m, precision);
  return out;
}

class SURidentical {
  public: 
    SURidentical(const mat&, const mat&, const mat&, 
                 const vec&, const mat&, int, int, int);
    double logML();
    mat g_draws, Omega_inv_draws;
  private:
    int r1, T, D, K, p, n_vech, j;
    vec G0_inv_g0, gbar, g;
    mat XX, XY, R0_inv, G0_inv;
    mat resid, RT, GT_inv, Omega_inv, RT_draws;
};
//Class constructor
SURidentical::SURidentical(const mat& X, const mat& Y,
                           const mat& G0, const vec& g0,
                           const mat& R0, int r0,
                           int n_draws = 1000,
                           int burn_in = 1000){
  T = Y.n_rows;
  D = Y.n_cols;
  K = X.n_cols;
  p = K * D;
  n_vech = (D + 1) * D / 2;
  
  Omega_inv_draws.zeros(n_vech, n_draws);
  g_draws.zeros(p, n_draws);
  RT_draws.zeros(n_vech, n_draws);
  XX = X.t() * X;
  XY = X.t() * Y;
  r1 = r0 + T; 
  R0_inv = inv_sympd(R0);
  G0_inv = inv_sympd(G0);
  G0_inv_g0 = solve(symmatu(G0), g0);
  
  Omega_inv = eye(D, D);
  
  for(int i = 0; i < (n_draws + burn_in); i++){
    
    GT_inv = G0_inv + kron(Omega_inv, XX);
    gbar = solve(GT_inv, G0_inv_g0 + vectorise(XY * Omega_inv)); 
    g = draw_normal(gbar, GT_inv);
    resid = Y - X * reshape(g, K, D);
    RT = inv_sympd(R0_inv, resid.t() * resid);
    Omega_inv = draw_wishart(r1, RT);
    
    if(i >= burn_in){
      j = i - burn_in;
      RT_draws.col(j) = vech(RT);
      Omega_inv_draws.col(j) = vech(Omega_inv);
      g_draws.col(j) = g;
    }
  }
}
//Member function to calculate marginal likelihood
double SURidentical::logML(){
  //calculate posterior means
  vec gstar = mean(g_draws, 1);
  mat Omega_inv_star = devech(mean(Omega_inv_draws, 1), D);
  //Still need to write body of this function...
  double prior_contrib, like_contrib, post_contrib;
  return prior_contrib + like_contrib + post_contrib;
}