// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::interfaces(r, cpp)]]

// ---------------------
// Miscellaneous
// ---------------------
arma::vec dif_exp(
    arma::vec& theta,
    double x
){
  arma::vec dx(2);
  dx(0) = std::exp(theta(0) + theta(1) * x);
  dx(1) = x * dx(0);
  
  return dx;
}

arma::vec dif_exp_1(
    arma::vec& theta,
    arma::vec& x
){
  return arma::exp(theta(0) + theta(1) * x);
}

arma::vec dif_exp_2(
    arma::vec& theta,
    arma::vec& x
){
  return x % dif_exp_1(theta,x);
}

arma::vec dif_tanh(
    arma::vec& theta,
    double x
){
  arma::vec dx(2);
  double t1;
  t1 = std::tanh(theta(0) + theta(1) * x);
  dx(0) = 0.1e1 - t1 * t1;
  dx(1) = x * dx(0);
  
  return dx;
}

double dif_tanh_1(
    double t
){
  return 0.1e1 - t * t;
}

arma::vec dif_tanh_1(
    arma::vec& theta,
    arma::vec& x
){
  unsigned int K(x.n_elem);
  arma::vec t1(K);
  
  for(unsigned int k(0); k<K; k++){
    t1(k) = std::tanh(theta(0) + theta(1) * x(k));
  }
  
  return 0.1e1 - t1 % t1;
}

double dif_tanh_2(
    double t,
    double x
){
  return x * dif_tanh_1(t);
}

arma::vec dif_tanh_2(
    arma::vec& theta,
    arma::vec& x
){
  return x % dif_tanh_1(theta,x);
}

// ---------------------
// Wavelet variance of the latent processes (data-dependent)
// ---------------------
// White Noise
arma::vec nu_wn(
    arma::vec& theta,
    double x,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  double gamma2;
  arma::vec nu(n);
  
  // Computation
  gamma2 = std::exp(theta(0) + theta(1) * x);
  nu = gamma2 / tau;
  
  return nu;
}

arma::mat nu_wn(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  unsigned int K(x.n_elem);
  arma::mat nu(n,K);
  arma::vec v1(n);
  
  // Computation
  v1 = std::exp(theta(0))  / tau;
  for(unsigned int j(0); j<K; j++){
    nu.col(j) = v1 * std::exp(theta(1) * x(j));
  }
  
  return nu;
}

arma::vec dif_nu_wn(
    arma::vec& tau
){
  return 0.1e1 / tau;
}

arma::mat jac_nu_wn_1(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  return dif_nu_wn(tau) * arma::trans(dif_exp_1(theta,x));
}

arma::mat jac_nu_wn_2(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  return dif_nu_wn(tau) * arma::trans(dif_exp_2(theta,x));
}

// Quantization noise
arma::vec nu_qn(
    arma::vec& theta,
    double x,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  double q2;
  arma::vec nu(n);
  
  // Computation
  q2 = 0.6e1 * std::exp(theta(0) + theta(1) * x);
  nu = q2 / tau / tau;
  
  return nu;
}

arma::mat nu_qn(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  unsigned int K(x.n_elem);
  arma::vec v1(n);
  arma::mat nu(n,K);
  
  // Computation
  v1 = 0.6e1 * std::exp(theta(0)) / tau / tau;
  for(unsigned int j(0); j<K; j++){
    nu.col(j) = v1 * std::exp(theta(1) * x(j));
  }
  
  return nu;
}

arma::vec dif_nu_qn(
    arma::vec& tau
){
  return 0.6e1 / tau / tau;
}

arma::mat jac_nu_qn_1(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  return dif_nu_qn(tau) * arma::trans(dif_exp_1(theta,x));
}

arma::mat jac_nu_qn_2(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  return dif_nu_qn(tau) * arma::trans(dif_exp_2(theta,x));
}

// Autoregressive process of order 1
arma::vec nu_ar1(
    arma::vec& theta,
    double x,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  double eta2;
  double phi;
  arma::vec phi_tau(n);
  arma::vec nu(n);
  double t1,t2;
  
  // Computation
  phi = std::tanh(theta(0) + theta(1) * x);
  eta2 = std::exp(theta(2) + theta(3) * x);
  
  for(unsigned int j = 0; j < n; j++){
    phi_tau(j) = std::pow(phi, tau(j));
  }
  
  t2 = phi - 0.1e1;
  t1 = eta2 / t2 / t2 / t2;
  t2 += 0.2e1;
  t1 /= t2;
  
  nu = t1 * (tau * phi * phi - tau + 0.6e1 * phi + 0.2e1 * phi * phi_tau - 0.8e1 * phi * arma::sqrt(phi_tau)) / tau / tau;
  
  return nu;
}

arma::mat nu_ar1(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  unsigned int K(x.n_elem);
  double eta2;
  double phi;
  arma::vec phi_tau(n);
  arma::vec phi_tau2(n);
  arma::mat nu(n,K);
  arma::vec v1(n);
  double t1,t2;
  
  // Computation
  for(unsigned int j(0); j<K; j++){
    phi = std::tanh(theta(0) + theta(1) * x(j));
    eta2 = std::exp(theta(2) + theta(3) * x(j));
    for(unsigned int k(0); k < n; k++){
      phi_tau(k) = std::pow(phi, tau(k));
    }
    phi_tau2 = arma::sqrt(phi_tau);
    v1 = tau % tau;
    t2 = phi - 0.1e1;
    t1 = eta2 / t2 / t2 / t2;
    t2 += 0.2e1;
    t1 /= t2;
    nu.col(j) = t1 * (tau * phi * phi - tau + 0.6e1 * phi + 0.2e1 * phi * phi_tau - 0.8e1 * phi * phi_tau2) / v1;
  }
  
  
  return nu;
}

arma::vec dif_nu_ar1_phi(
    arma::vec& lambda,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  double eta2;
  double phi;
  arma::vec phi_tau(n);
  arma::vec phi_tau2(n);
  arma::vec dn(n);
  double t1,t2,t3;
  arma::vec v1(n);
  arma::vec v2(n);
  arma::vec v3(n);
  
  // Computation
  phi = lambda(0);
  eta2 = lambda(1);
  
  t2 = phi - 0.1e1;
  t3 = phi + 0.1e1;
  t1 = 0.2e1 * eta2 / t2 / t2 / t2 / t2 / t3 / t3;
  t2 = phi * phi; 
  
  for(unsigned int j = 0; j < n; j++){
    phi_tau(j) = std::pow(phi, tau(j));
  }
  phi_tau2 = arma::sqrt(phi_tau);
  v1 = tau * t2 - tau;
  v2 = tau % tau;
  v3 = phi_tau - 0.4e1 * phi_tau2 + 0.3e1;
  phi_tau -= 0.2e1 * phi_tau2 + t3;
  t3 = 0.3e1 * t2 + 0.2e1 * phi + 0.1e1;
  
  dn = t1 * v1 % phi_tau / v2 - t1 * t3 * v3 / v2;
  
  return dn;
}

arma::mat jac_nu_ar1_phi_1(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int K(x.n_elem);
  unsigned int n(tau.n_elem);
  arma::mat jac(n,K);
  arma::vec par(2);
  
  for(unsigned int k(0); k<K; k++){
    par(0) = std::tanh(theta(0) + theta(1) * x(k));
    par(1) = std::exp(theta(2) + theta(3) * x(k));
    jac.col(k) = dif_nu_ar1_phi(par,tau) * dif_tanh_1(par(0));
  }
  
  return jac;
}

arma::mat jac_nu_ar1_phi_2(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int K(x.n_elem);
  unsigned int n(tau.n_elem);
  arma::mat jac(n,K);
  arma::vec par(2);
  
  for(unsigned int k(0); k<K; k++){
    par(0) = std::tanh(theta(0) + theta(1) * x(k));
    par(1) = std::exp(theta(2) + theta(3) * x(k));
    jac.col(k) = dif_nu_ar1_phi(par,tau) * dif_tanh_2(par(0),x(k));
  }
  
  return jac;
}

arma::vec dif_nu_ar1_eta2(
    arma::vec& lambda,
    arma::vec& tau
){
  unsigned int n(tau.n_elem);
  // double eta2;
  double phi;
  arma::vec phi_tau(n);
  arma::vec phi_tau2(n);
  arma::vec dn(n);
  double t1,t2,t3;
  arma::vec v1(n);
  arma::vec v2(n);
  arma::vec v3(n);
  
  // Computation
  phi = lambda(0);
  // eta2 = lambda(1);
  
  for(unsigned int j = 0; j < n; j++){
    phi_tau(j) = std::pow(phi, tau(j));
  }
  phi_tau2 = arma::sqrt(phi_tau);
  v1 = tau * phi * phi - tau;
  v2 = tau % tau;
  v3 = phi_tau - 0.4e1 * phi_tau2 + 0.3e1;
  
  t2 = phi - 0.1e1;
  t3 = phi + 0.1e1;
  t1 = 0.1e1 / t2 / t2 / t2 / t3;
  
  dn = t1 * v1 / v2 + t1 * 0.2e1 * phi * v3 / v2;
  
  return dn;
}

arma::mat jac_nu_ar1_eta2_1(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int K(x.n_elem);
  unsigned int n(tau.n_elem);
  arma::mat jac(n,K);
  arma::vec par(2);
  
  for(unsigned int k(0); k<K; k++){
    par(0) = std::tanh(theta(0) + theta(1) * x(k));
    par(1) = std::exp(theta(2) + theta(3) * x(k));
    jac.col(k) = dif_nu_ar1_eta2(par,tau) * par(1);
  }
  
  return jac;
}

arma::mat jac_nu_ar1_eta2_2(
    arma::vec& theta,
    arma::vec& x,
    arma::vec& tau
){
  unsigned int K(x.n_elem);
  unsigned int n(tau.n_elem);
  arma::mat jac(n,K);
  arma::vec par(2);
  
  for(unsigned int k(0); k<K; k++){
    par(0) = std::tanh(theta(0) + theta(1) * x(k));
    par(1) = std::exp(theta(2) + theta(3) * x(k));
    jac.col(k) = dif_nu_ar1_eta2(par,tau) * par(1) * x(k);
  }
  
  return jac;
}

// ---------------------
// Objective function dynamic GMWM estimator
// ---------------------
//' Objective function for estimating dynamic GMWM parameters
//' 
//' @param theta  A vector of parameter
//' @param nu_hat A matrix of wavelet variances
//' @param x      A vector of data
//' @param tau    A vector containing the scales
//' @param Omega  A matrix of weights
//' @param WN     A boolean inidicating whether a gaussian White Noise is 
//' included in the model
//' @param QN     A boolean inidicating whether a Quantization Noise is included 
//' in the model
//' @param AR1    A non-negative integer indicating the number of Autoregressive 
//' process(es) of order 1 to be included in the model
//' @return Returns a double \deqn{\frac{1}{2}\frac{1}{K}\sum_{k=1}^K\lVert\hat{\nu}_k - \nu(\theta,x_k)\rVert^2_{\Omega}}
//' @export
// [[Rcpp::export]]
double of_dyn(
    arma::vec& theta,
    arma::mat& nu_hat,
    arma::vec& x,
    arma::vec& tau,
    arma::mat& Omega,
    bool WN,
    bool QN,
    unsigned int AR1
){
  unsigned int K(x.n_elem); // number of bins
  unsigned int n(tau.n_elem); // number of scale
  arma::mat nu(n,K,arma::fill::zeros);
  double of(0.0);
  unsigned int d(0);
  
  if (n != nu_hat.n_rows) {
    Rcpp::stop("Mismatch: nu_hat should have the same number of rows as the size of tau!");
  }
  
  if (K != nu_hat.n_cols) {
    Rcpp::stop("Mismatch: nu_hat should have the same number of columns as the size of x!");
  }
  
  if ((n != Omega.n_cols) || (n != Omega.n_rows)) {
    Rcpp::stop("The matrix Omega has not the right dimension");
  }
  
  if (WN) {
    d += 2;
  }
  
  if (QN) {
    d += 2;
  }
  
  if (AR1 != 0) {
    d += AR1 * 4;
  }
  
  if (d != theta.n_elem) {
    Rcpp::stop("You have not provided the right number of parameter theta.");
  }
  
  d = 0;
  
  if (WN) {
    arma::vec t1(2);
    t1 = theta.subvec(d,d+1);
    nu += nu_wn(t1,x,tau);
    d += 2;
  }
  
  if (QN) {
    arma::vec t1(2);
    t1 = theta.subvec(d,d+1);
    nu += nu_qn(t1,x,tau);
    d += 2;
  }
  
  if (AR1 != 0) {
    arma::vec t1(4);
    if (AR1 == 1) {
      t1 = theta.subvec(d,d+3);
      nu += nu_ar1(t1,x,tau);
    } else {
      for(unsigned int k(0); k < AR1; k++){
        t1 = theta.subvec(d,d+3);
        nu += nu_ar1(t1,x,tau);
        d += 4;
      }
    }
  }
  
  nu -= nu_hat;
  for(unsigned int k(0); k<K; k++){
    of += arma::dot(nu.col(k), Omega * nu.col(k));
  }
  
  return of / 0.2e1 / K;
}

//' Gradient of the objective function for estimating dynamic GMWM parameters
//' 
//' @param theta  A vector of parameter
//' @param nu_hat A matrix of wavelet variances
//' @param x      A vector of data
//' @param tau    A vector containing the scales
//' @param Omega  A matrix of weights
//' @param WN     A boolean inidicating whether a gaussian White Noise is 
//' included in the model
//' @param QN     A boolean inidicating whether a Quantization Noise is included 
//' in the model
//' @param AR1    A non-negative integer indicating the number of Autoregressive 
//' process(es) of order 1 to be included in the model
//' @return Returns a vector \deqn{\frac{1}{K}\sum_{k=1}^K(\hat{\nu}_k - \nu(\theta,x_k))\Omega\frac{\partial}{\partial\theta}\nu(\theta,x_k)}
//' @export
// [[Rcpp::export]]
arma::vec grad_of_dyn(
    arma::vec& theta,
    arma::mat& nu_hat,
    arma::vec& x,
    arma::vec& tau,
    arma::mat& Omega,
    bool WN,
    bool QN,
    unsigned int AR1
){
  unsigned int K(x.n_elem); // number of bins
  unsigned int n(tau.n_elem); // number of scale
  unsigned int p(theta.n_elem); // number of parameters
  arma::mat nu(n,K,arma::fill::zeros);
  arma::mat m(n,K);
  arma::cube d_nu(n,K,p);
  arma::vec grad(p);
  unsigned int d(0);

  if (n != nu_hat.n_rows) {
    Rcpp::stop("Mismatch: nu_hat should have the same number of rows as the size of tau!");
  }

  if (K != nu_hat.n_cols) {
    Rcpp::stop("Mismatch: nu_hat should have the same number of columns as the size of x!");
  }

  if ((n != Omega.n_cols) || (n != Omega.n_rows)) {
    Rcpp::stop("The matrix Omega has not the right dimension");
  }

  if (WN) {
    d += 2;
  }

  if (QN) {
    d += 2;
  }

  if (AR1 != 0) {
    d += AR1 * 4;
  }

  if (d != p) {
    Rcpp::stop("You have not provided the right number of parameter theta.");
  }

  d = 0;

  if (WN) {
    arma::vec t1(2);
    t1 = theta.subvec(d,d+1);
    nu += nu_wn(t1,x,tau);
    d_nu.slice(d) = jac_nu_wn_1(t1,x,tau);
    d_nu.slice(d+1) = jac_nu_wn_2(t1,x,tau);
    d += 2;
  }

  if (QN) {
    arma::vec t1(2);
    t1 = theta.subvec(d,d+1);
    nu += nu_qn(t1,x,tau);
    d_nu.slice(d) = jac_nu_qn_1(t1,x,tau);
    d_nu.slice(d+1) = jac_nu_qn_2(t1,x,tau);
    d += 2;
  }

  if (AR1 != 0) {
    arma::vec t1(4);
    if( AR1 == 1) {
      t1 = theta.subvec(d,d+3);
      nu += nu_ar1(t1,x,tau);
      d_nu.slice(d) = jac_nu_ar1_phi_1(t1,x,tau);
      d_nu.slice(d+1) = jac_nu_ar1_phi_2(t1,x,tau);
      d_nu.slice(d+2) = jac_nu_ar1_eta2_1(t1,x,tau);
      d_nu.slice(d+3) = jac_nu_ar1_eta2_2(t1,x,tau);
    } else {
      for(unsigned int k(0); k < AR1; k++){
        t1 = theta.subvec(d,d+3);
        nu += nu_ar1(t1,x,tau);
        d_nu.slice(d) = jac_nu_ar1_phi_1(t1,x,tau);
        d_nu.slice(d+1) = jac_nu_ar1_phi_2(t1,x,tau);
        d_nu.slice(d+2) = jac_nu_ar1_eta2_1(t1,x,tau);
        d_nu.slice(d+3) = jac_nu_ar1_eta2_2(t1,x,tau);
        d += 4;
      }
    }
  }

  nu -= nu_hat;
  for(unsigned int j(0); j<p; j++){
    m = d_nu.slice(j);
    for(unsigned int k(0); k<K; k++){
      grad(j) += arma::dot(m.col(k), Omega * nu.col(k));
    }
  }

  return grad / K;
}
