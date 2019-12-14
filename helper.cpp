#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;



// Functions for other helper functions
// ####################################

double dbetabinom(int n, int k, double alpha, double beta) {
  return(exp(Rf_lchoose(n,k) + Rf_lbeta(k+alpha,n-k+beta) - Rf_lbeta(alpha,beta)));
}

double expit(double x) {
  return(1/(1+exp(-x)));
}





// Functions for linear model
// ##########################

// [[Rcpp::export]]
NumericMatrix est_params(NumericMatrix clpost, NumericMatrix friends_x, NumericVector y) {
  // Estimate simple linear regression parameters assuming a distribution over covariate x for each exposure condition
  NumericMatrix params(3,4);

  int xlen = friends_x.ncol();

  int i, j;
  double sum_x, sum_y, sum_xx, sum_xy, sum_resid, denom;
  NumericVector col_c, col_x, resid;

  for(i = 0; i < 4; i++) {
    sum_x = 0;
    sum_y = 0;
    sum_xx = 0;
    sum_xy = 0;
    sum_resid = 0;
    denom = 0;

    for(j = 0; j < xlen; j++) {
      col_c = clpost(_,j + i*xlen);
      col_x = friends_x(_,j);
      sum_x += sum(col_x*col_c);
      sum_y += sum(y*col_c);
      sum_xx += sum(col_x*col_x*col_c);
      sum_xy += sum(col_x*y*col_c);
      denom += sum(col_c);
    }

    // slope
    params(0,i) = (sum_xy - sum_x*sum_y/denom)/(sum_xx - sum_x*sum_x/denom);
    // intercept
    params(1,i) = sum_y/denom - params(0,i)*sum_x/denom;

    for(j = 0; j < xlen; j++) {
      resid = y - params(1,i) - params(0,i)*friends_x(_,j);
      sum_resid += sum(clpost(_,j + i*xlen) * resid * resid);
    }

    // sigma
    params(2,i) = sqrt(sum_resid/denom);
  }

  return(params);
}

// [[Rcpp::export]]
NumericMatrix est_params_eqvar(NumericMatrix clpost, NumericMatrix friends_x, NumericVector y) {
    // Estimate simple linear regression parameters assuming a distribution over covariate x for each exposure condition
    // assuming the models for each exposure condition have equal variance
  NumericMatrix params(3,4);

  int xlen = friends_x.ncol();

  int i, j;
  double sum_x, sum_y, sum_xx, sum_xy, sum_resid, denom, sum_denom;
  NumericVector col_c, col_x, resid;

  sum_resid = 0;
  sum_denom = 0;

  for(i = 0; i < 4; i++) {
    sum_x = 0;
    sum_y = 0;
    sum_xx = 0;
    sum_xy = 0;
    denom = 0;

    for(j = 0; j < xlen; j++) {
      col_c = clpost(_,j + i*xlen);
      col_x = friends_x(_,j);
      sum_x += sum(col_x*col_c);
      sum_y += sum(y*col_c);
      sum_xx += sum(col_x*col_x*col_c);
      sum_xy += sum(col_x*y*col_c);
      denom += sum(col_c);
    }

    // slope
    params(0,i) = (sum_xy - sum_x*sum_y/denom)/(sum_xx - sum_x*sum_x/denom);
    // intercept
    params(1,i) = sum_y/denom - params(0,i)*sum_x/denom;

    for(j = 0; j < xlen; j++) {
      resid = y - params(1,i) - params(0,i)*friends_x(_,j);
      sum_resid += sum(clpost(_,j + i*xlen) * resid * resid);
    }
    sum_denom += denom;
  }

  // sigma
  for(i = 0; i < 4; i++) {
    params(2,i) = sqrt(sum_resid/sum_denom);
  }

  return(params);
}



// [[Rcpp::export]]
NumericMatrix est_clprior(int min_missing, int max_missing, NumericVector itt, NumericVector itt2, int N2, NumericVector num_treated,
                          NumericVector trt_friends, NumericVector num_friends, double p_prob, double q_prob, double alpha, double beta) {
  // Estimate prior of true exposure condition and degree distribution
  int xlen = max_missing - min_missing + 1;
  int ylen = itt.size();
  int Nx = N2 - 1;

  NumericMatrix clprior(4*xlen, ylen);
  //Distribution of treated and non-treated friends
  double prob_friend_t[(int) max(num_friends) + max_missing];
  double prob_friend_nt[(int) max(num_friends) + max_missing];

  int i,j,k,num_dropped;
  int num_other_trt, delta_friend, max_friend;
  double temp_sum, temp_sum2;

  for(i = 0; i < ylen; i++) {
    max_friend = num_friends[i] + max_missing;
    num_other_trt = num_treated[i] - itt2[i];

    for(j = 0; j <= max_friend; j++) {
      temp_sum = 0;
      temp_sum2 = 0;
      for(num_dropped = 0; num_dropped <= j; num_dropped++) {
        if(num_other_trt >= j) {
          temp_sum += dbetabinom(num_other_trt,j,alpha,beta) * Rf_dbinom(num_dropped, j, p_prob, 0) *
                      Rf_dbinom(trt_friends[i] - j + num_dropped, num_other_trt - j, q_prob, 0);
        }
        temp_sum2 += dbetabinom(Nx - num_other_trt,j,alpha,beta) * Rf_dbinom(num_dropped, j, p_prob, 0) *
                     Rf_dbinom(num_friends[i] - trt_friends[i] - j + num_dropped, Nx - num_other_trt - j, q_prob, 0);
      }

      // Probability of actually having j treated / non-treated friends,
      // calculated by marginalizing over the number of friends that are not observed
      prob_friend_t[j] = temp_sum;
      prob_friend_nt[j] = temp_sum2;
    }

    for(j = 0; j < max_friend; j++) {
      for(k = 0; k < max_friend; k++) {
        // Difference between the number of actual friends and the number of observed friends
        delta_friend = j + k - num_friends[i];

        if((delta_friend + num_friends[i] >= 0) && (delta_friend >= min_missing) && (delta_friend <= max_missing)) {
          // Lexographical ordering of exposure condition and number of friends
          clprior(delta_friend - min_missing + (2*itt[i] + (j>0))*xlen,i) += prob_friend_t[j]*prob_friend_nt[k];
        }
      }
    }
  }

  return(clprior);
}



// [[Rcpp::export]]
NumericMatrix est_clpost(NumericMatrix clprior, NumericMatrix friends_x, NumericVector y, NumericMatrix params) {
  // Estimate posterior of true exposure condition and degree distribution assuming a linear model
  NumericMatrix clpost(clprior.nrow(), clprior.ncol());

  int xlen = friends_x.ncol();

  int i,j;
  for(i = 0; i < 4; i++) {
    for(j = 0; j < xlen; j++) {
       clpost(_,i*xlen + j) = clprior(_,i*xlen + j) * dnorm(y - params(1,i) - params(0,i)*friends_x(_,j), 0.0, params(2,i));
    }
  }

  return(clpost);
}



// [[Rcpp::export]]
double corr_objective(NumericVector x, NumericMatrix post_mat, NumericVector num_friends, int N2, int min_missing, double alpha, double beta) {
  // Objective function for estimating p and q (x)
  double obj = 0;

  int ylen = post_mat.ncol();
  int xlen = post_mat.nrow();

  int i,j,k,num_friends_true;
  NumericVector temp_sum(xlen);
  double cum_sum;

  // Iterate over subjects i
  for(i = 0; i < ylen; i++) {
    // Iterate over number of friends
    for(j = 0; j < xlen; j++) {
      temp_sum[j] = 0;
      if(post_mat(j,i)) {
        num_friends_true = num_friends[i] + min_missing + j;
        for(k = std::max(0,j+min_missing); k <= num_friends_true; k++) {
          temp_sum[j] += dbetabinom(N2-1,num_friends_true,alpha,beta) *
                         Rf_choose(num_friends_true, k) * pow(x[0], k) * pow(1-x[0],num_friends_true-k) *
                         Rf_choose(N2 - num_friends_true - 1, k - j - min_missing) *
                         pow(x[1], k - j - min_missing) * pow(1-x[1], N2 - num_friends_true - 1 + j + min_missing - k);
        }
      }
    }

    cum_sum = sum(temp_sum);
    for(j = 0; j < xlen; j++) {
      if(post_mat(j,i)) {
        obj += post_mat(j,i)*(log(temp_sum[j]) - log(cum_sum));
      }
    }
  }

  return(-obj);
}

// [[Rcpp::export]]
NumericVector corr_gradient(NumericVector x, NumericMatrix post_mat, NumericVector num_friends, int N2, int min_missing, double alpha, double beta) {
  // Gradient for objective function corr_objective
  double obj1 = 0;
  double obj2 = 0;

  int ylen = post_mat.ncol();
  int xlen = post_mat.nrow();

  int i,j,k, num_friends_true;
  double cum_numer1, cum_numer2, cum_denom;
  NumericVector numer1(xlen);
  NumericVector numer2(xlen);
  NumericVector denom(xlen);

  // Iterate over subjects i
  for(i = 0; i < ylen; i++) {
    // Iterate over number of friends
    for(j = 0; j < xlen; j++) {
      numer1[j] = 0;
      numer2[j] = 0;
      denom[j] = 0;
      if(post_mat(j,i)) {
        num_friends_true = num_friends[i] + min_missing + j;
        for(k = std::max(0,j+min_missing); k <= num_friends_true; k++) {
          numer1[j] += (k-x[0]*(num_friends_true)) * dbetabinom(N2-1,num_friends_true,alpha,beta) *
                       Rf_choose(num_friends_true, k) * pow(x[0], k-1) * pow(1-x[0],num_friends_true-k-1) *
                       Rf_choose(N2 - num_friends_true - 1, k - j - min_missing) *
                      pow(x[1], k - j - min_missing) * pow(1-x[1],N2 - num_friends_true - 1 + j + min_missing - k);
          numer2[j] += (k - j - min_missing - (N2 - num_friends_true - 1)*x[1]) * dbetabinom(N2-1,num_friends_true,alpha,beta) *
                       Rf_choose(num_friends_true, k) * pow(x[0], k) * pow(1-x[0],num_friends_true-k) *
                       Rf_choose(N2 - num_friends_true - 1, k - j - min_missing) *
                       pow(x[1], k - j - min_missing - 1) * pow(1-x[1],N2 - num_friends_true - 1 + j + min_missing - k - 1);
          denom[j] += dbetabinom(N2-1,num_friends_true,alpha,beta) *
                      Rf_choose(num_friends_true, k) * pow(x[0], k) * pow(1-x[0],num_friends_true-k) *
                      Rf_choose(N2 - num_friends_true - 1, k - j - min_missing) *
                      pow(x[1], k - j - min_missing) * pow(1-x[1],N2 - num_friends_true - 1 + j + min_missing - k);
        }
      }
    }

    cum_numer1 = sum(numer1);
    cum_numer2 = sum(numer2);
    cum_denom = sum(denom);
    for(j = 0; j < xlen; j++) {
      if(post_mat(j,i)) {
        obj1 += post_mat(j,i)*(numer1[j]/denom[j] - cum_numer1/cum_denom);
        obj2 += post_mat(j,i)*(numer2[j]/denom[j] - cum_numer2/cum_denom);
      }
    }
  }

  NumericVector obj(2);
  obj[0] = -obj1;
  obj[1] = -obj2;

  return(obj);
}










// Functions for logistic model
// ############################

// [[Rcpp::export]]
double param_binom_objective(NumericVector x, NumericMatrix post_mat, NumericMatrix num_friends_x, NumericVector n, NumericVector a) {
  // Objective function for logistic likelihood to be maximized w.r.t. parameters x
  double obj = 0;

  int N = post_mat.nrow();
  int xlen = post_mat.ncol();

  int i,j;

  for(i = 0; i < N; i++) {
    for(j = 0; j < xlen; j++) {
      if(post_mat(i,j)) {
        obj -= post_mat(i,j)*R::dbinom(a[i],n[i],expit(x[0] + x[1]*num_friends_x(i,j)),1);
      }
    }
  }

  return(obj);
}



// [[Rcpp::export]]
NumericMatrix est_clpost_binom(NumericMatrix clprior, NumericMatrix friends_x, NumericVector n, NumericVector a, NumericMatrix params) {
  // Estimate posterior of true exposure condition and degree distribution assuming a logistic model
  NumericMatrix clpost(clprior.nrow(), clprior.ncol());

  int xlen = friends_x.ncol();
  int N = friends_x.nrow();

  int i,j,k;
  for(i = 0; i < 4; i++) {
    for(j = 0; j < xlen; j++) {
      for(k = 0; k < N; k++) {
        clpost(k,i*xlen + j) = clprior(k,i*xlen + j) * R::dbinom(a[k],n[k],expit(params(0,i) + params(1,i)*friends_x(k,j)),0);
      }
    }
  }

  return(clpost);
}










// Functions for the logistic model including in- and out-of-village covariates
// ############################################################################

// [[Rcpp::export]]
NumericMatrix est_clprior_io(int min_missing, int max_missing, NumericVector itt, NumericVector itt2, NumericVector N_in, NumericVector N_out,
                             NumericVector num_treated_in, NumericVector num_treated_out, NumericVector trt_friends_in, NumericVector trt_friends_out,
                            NumericVector num_friends_in, NumericVector num_friends_out,
                            double p_prob_in, double q_prob_in, double p_prob_out, double q_prob_out,
                            double alpha_in, double beta_in, double alpha_out, double beta_out) {
  // Estimate prior of true exposure condition and degree distribution
  int xlen = max_missing - min_missing + 1;
  int ylen = itt.size();

  NumericMatrix clprior(4*xlen, ylen);

  //Distribution of treated and non-treated friends (in-village, out-of-village, and combined)
  double prob_friend_in_t[(int) max(num_friends_in) + max_missing];
  double prob_friend_in_nt[(int) max(num_friends_in) + max_missing];
  double prob_friend_out_t[(int) max(num_friends_out) + max_missing/2];
  double prob_friend_out_nt[(int) max(num_friends_out) + max_missing/2];
  double prob_friend_t[(int) max(num_friends_in + num_friends_out) + max_missing];
  double prob_friend_nt[(int) max(num_friends_in + num_friends_out) + max_missing];

  int i,j,k,num_dropped;
  int Nx, num_other_trt, delta_friend, max_friend, max_friend2, max_friend3;
  double temp_sum, temp_sum2;

  for(i = 0; i < ylen; i++) {
    Nx = N_in[i];
    max_friend = num_friends_in[i] + max_missing > Nx ? Nx : num_friends_in[i] + max_missing;
    num_other_trt = num_treated_in[i];

    for(j = 0; j <= max_friend; j++) {
      temp_sum = 0;
      temp_sum2 = 0;
      for(num_dropped = 0; num_dropped <= j; num_dropped++) {
        if(num_other_trt >= j) {
          temp_sum += exp(Rf_dbinom(num_dropped, j, p_prob_in, 1) + Rf_dbinom(trt_friends_in[i] - j + num_dropped, num_other_trt - j, q_prob_in, 1) -
                          Rf_dbinom(trt_friends_in[i], num_other_trt, q_prob_in, 1));
        }
        if(Nx - num_other_trt >= j) {
          temp_sum2 += exp(Rf_dbinom(num_dropped, j, p_prob_in, 1) +
                           Rf_dbinom(num_friends_in[i] - trt_friends_in[i] - j + num_dropped, Nx - num_other_trt - j, q_prob_in, 1) -
                           Rf_dbinom(num_friends_in[i] - trt_friends_in[i], Nx - num_other_trt, q_prob_in, 1));
        }
      }

      // Probability of actually having j treated / non-treated friends in village,
      // calculated by marginalizing over the number of friends that are not observed
      prob_friend_in_t[j] = temp_sum * dbetabinom(num_other_trt,j,alpha_in,beta_in);
      prob_friend_in_nt[j] = temp_sum2 * dbetabinom(Nx - num_other_trt,j,alpha_in,beta_in);
    }

    Nx = N_out[i];
    max_friend2 = num_friends_out[i] + max_missing/2;
    num_other_trt = num_treated_out[i];

    for(j = 0; j <= max_friend2; j++) {
      temp_sum = 0;
      temp_sum2 = 0;
      for(num_dropped = 0; num_dropped <= j; num_dropped++) {
        if(num_other_trt >= j) {
          temp_sum += Rf_dbinom(num_dropped, j, p_prob_out, 0) * Rf_dbinom(trt_friends_out[i] - j + num_dropped, num_other_trt - j, q_prob_out, 0);
        }
        if(Nx - num_other_trt >= j) {
          temp_sum2 += Rf_dbinom(num_dropped, j, p_prob_out, 0) *
                       Rf_dbinom(num_friends_out[i] - trt_friends_out[i] - j + num_dropped, Nx - num_other_trt - j, q_prob_out, 0);
        }
      }

      // Probability of actually having j treated / non-treated friends out of  village,
      // calculated by marginalizing over the number of friends that are not observed
      prob_friend_out_t[j] = temp_sum * dbetabinom(num_other_trt,j,alpha_out,beta_out);
      prob_friend_out_nt[j] = temp_sum2 * dbetabinom(Nx - num_other_trt,j,alpha_out,beta_out);
    }

    max_friend3 = num_friends_in[i] + num_friends_out[i] + max_missing;
    for(j = 0; j < max_friend3; j++) {
      prob_friend_t[j] = 0;
      prob_friend_nt[j] = 0;
    }

    for(j = 0; j <= max_friend; j++) {
      for(k = 0; k <= max_friend2; k++) {
        if(j + k < max_friend3) {
          // Probability of j+k treated / non-treated friends
          prob_friend_t[j+k] += prob_friend_in_t[j]*prob_friend_out_t[k];
          prob_friend_nt[j+k] += prob_friend_in_nt[j]*prob_friend_out_nt[k];
        }
      }
    }

    for(j = 0; j < max_friend3; j++) {
      for(k = 0; k < max_friend3; k++) {
        delta_friend = j + k - num_friends_in[i] - num_friends_out[i];
        if((delta_friend >= min_missing) && (delta_friend <= max_missing)) {
          // Lexographical ordering of exposure condition and number of friends
          clprior(delta_friend - min_missing + (2*itt[i] + (j>0))*xlen,i) += prob_friend_t[j]*prob_friend_nt[k];
        }
      }
    }
  }

  return(clprior);
}



// [[Rcpp::export]]
double corr_objective_io(NumericVector x, NumericMatrix post_mat, NumericVector num_friends_in, NumericVector num_friends_out,
                       NumericVector N_in, NumericVector N_out, int min_missing, int max_missing,
                       double alpha_in, double beta_in, double alpha_out, double beta_out) {
  // Objective function for estimating p_in, q_in, p_out, and q_out (x)
  double obj = 0;

  int ylen = post_mat.ncol();
  int xlen = post_mat.nrow();

  int i,j,k,l,max_friend,max_friend2,Nx;
  NumericVector temp_sum(xlen);
  double cum_sum, temp_sum2;

  double prob_friend_in[(int) max(num_friends_in) + max_missing];
  double prob_friend_out[(int) max(num_friends_out) + max_missing];

  //Iterate over subjects j
  for(i = 0; i < ylen; i++) {
    for(j = 0; j < xlen; j++) {
      temp_sum[j] = 0;
    }

    Nx = N_in[i] - 1;
    max_friend = num_friends_in[i] + max_missing > Nx ? Nx : num_friends_in[i] + max_missing;

    for(j = 0; j <= max_friend; j++) {
      temp_sum2 = 0;
      for(k = 0; k <= j; k++) {
        temp_sum2 += exp(Rf_dbinom(k, j, x[0], 1) + Rf_dbinom(num_friends_in[i] - j + k, Nx - j, x[1], 1));
      }

      // Probability of j friends in-village, calculated by iterating over number of unobserved friends k
      prob_friend_in[j] = temp_sum2 * dbetabinom(Nx,j,alpha_in,beta_in);
    }

    Nx = N_out[i];
    max_friend2 = num_friends_out[i] + max_missing;

    for(j = 0; j <= max_friend2; j++) {
      temp_sum2 = 0;
      for(k = 0; k <= j; k++) {
        temp_sum2 += exp(Rf_dbinom(k, j, x[2], 1) + Rf_dbinom(num_friends_out[i] - j + k, Nx - j, x[3], 1));
      }

      // Probability of j friends in-village, calculated by iterating over number of unobserved friends k
      prob_friend_out[j] = temp_sum2 * dbetabinom(Nx,j,alpha_out,beta_out);
    }

    for(k = 0; k <= max_friend; k++) {
      for(l = 0; l <= max_friend2; l++) {
        if((k + l - num_friends_in[i] - num_friends_out[i] >= min_missing) && (k + l - num_friends_in[i] - num_friends_out[i] <= max_missing) &&
           (post_mat(k + l - num_friends_in[i] - num_friends_out[i] - min_missing,i) > 0)
        )
          temp_sum[k + l - num_friends_in[i] - num_friends_out[i] - min_missing] += prob_friend_in[k] * prob_friend_out[l];
      }
    }


    cum_sum = sum(temp_sum);
    for(j = 0; j < xlen; j++) {
      if((post_mat(j,i) > 0) && cum_sum > 0) {
        obj += post_mat(j,i)*(log(temp_sum[j]) - log(cum_sum));
      }
    }
  }

  return(-obj);
}
