// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "distributions.h"

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
colvec testy(colvec m){
  int p = m.n_elem;
  mat precision = eye(p, p);
  colvec out = draw_normal(m, precision);
  return out;
}