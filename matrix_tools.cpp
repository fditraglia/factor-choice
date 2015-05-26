// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdexcept>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
colvec vech(mat A){
  int m = A.n_rows;
  int n = A.n_cols;
  if(m != n){
    throw std::invalid_argument("received non-square matrix");
  }
  colvec out((n + 1) * n / 2, fill::zeros);
  int out_index = 0;
  for(int col = 0; col < n; col++){
    for(int row = col; row < n; row++){
      out(out_index) = A(row, col);
      out_index++;
    }
  }
  return out;
}


// [[Rcpp::export]]
mat devech(colvec v, int dim){
  if(v.n_elem != ((dim + 1) * dim / 2)){
    throw std::invalid_argument("dim and length v disagree");
  }
  mat out(dim, dim, fill::zeros);
  int v_index = 0;
  for(int col = 0; col < dim; col++){
    for(int row = col; row < dim; row++){
      out(row, col) = v(v_index);
      v_index++;
    }
  }
  return symmatl(out);
}

/*** R
M <- matrix(c(11, 12, 13,
              12, 22, 23,
              13, 23, 33), 3, 3, byrow = TRUE)
vech(M)
v <- drop(vech(M))
all.equal(devech(v, 3), M)
*/