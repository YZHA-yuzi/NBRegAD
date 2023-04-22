#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat get_predvals(arma::mat parm_post,
                       arma::mat X_all,
                       int ntotal){
  int npost = parm_post.n_rows;
  int nparm = parm_post.n_cols;
  arma::mat out;
  out.reshape(ntotal, npost);
  for(int i=0; i<npost; ++i){
    arma::mat parm_i = parm_post.row(i);
    parm_i.reshape(nparm, 1);
    arma::mat out1 = X_all*parm_i;
    out.col(i) = out1;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat get_count(arma::vec offset,
                    arma::mat pred_mat,
                    arma::vec overdisp_post){
  int ntotal = pred_mat.n_rows;
  for(int i=0; i<ntotal; ++i){
    double offset_i = offset(i);
    pred_mat.row(i) += offset_i;
  }
  arma::mat out1 = exp(pred_mat);
  for(int j=0; j<out1.n_cols; ++j){
    double xi_j = overdisp_post(j);
    out1.col(j) *= xi_j;
  }
  return out1;
}

