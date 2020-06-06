// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::export]]
Rcpp::List fastLm(const arma::colvec& y, const arma::mat& X) {
  int n = X.n_rows, k = X.n_cols;
  arma::colvec coef = solve(X, y);
  arma::colvec resid = y - X*coef;
  double sig2 = as_scalar(trans(resid) * resid/(n-k));
  arma::colvec stderrest = sqrt(sig2 * diagvec( inv(trans(X)*X)) );
  return Rcpp::List::create(Rcpp::_["coefficients"] = coef,
                            Rcpp::_["stderr"]       = stderrest,
                            Rcpp::_["df.residual"]  = n - k );
}