// load Rcpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;       // shorthand

// [[Rcpp::export]]
NumericMatrix bfmu_c(NumericMatrix x, int ite, double a = 2.1, double b = 1.1, double gamma_a = 1,
    double gamma_b = 1, double omega = 10, double omega_1 = 0.01) {

    int m = x.nrow();
    int n = x.ncol();

    arma::colvec alpha = rep(0.1, m);
    arma::rowvec sigma2 = rep(1.0, m);
    arma::rowvec z = rep(0.0, m);
    arma::rowvec p_star = rep(0.5, m);
    arma::rowvec lambda = rnorm(n, 0, 1);

    arma::rowvec debug = rep(0.0, m);
    arma::mat debug_matrix(ite, m);

    arma::mat dinv(m, m);
    dinv.fill(0.0);

    arma::mat alpha_matrix(ite, m);
    arma::mat lambda_matrix(ite, n);
    arma::mat sigma_matrix(ite, m);
    arma::mat z_matrix(ite, m);
    arma::mat p_matrix(ite, m);

    arma::mat x_arma(x.begin(), m, n, false);

    for (int k = 0; k < ite; k++) {
        for (int i = 0; i < m; i++) {
            double A = a + n/2;
            double B = b + 0.5*(arma::as_scalar(x_arma.row(i)*x_arma.row(i).t()) -
                2*arma::as_scalar(lambda*x_arma.row(i).t())*alpha[i] + 
                alpha[i]*arma::as_scalar(lambda*lambda.t())*alpha[i]);
            sigma2[i] = 1/R::rgamma(A, B);

            if(R::runif(0,1) < p_star[i]) {
                z[i] = 1;
            } else {
                z[i] = 0;
            }

            double gamma_1 = gamma_a + z[i];
            double gamma_2 = gamma_b + 1 - z[i];
            double p = R::rbeta(gamma_1, gamma_2);

            double lambda_ss = 0.0;
            for(int j = 0; j<n; j++) {
                lambda_ss += pow(lambda[j], 2);
            }
            double v_alpha_0 = 1/(1/omega + lambda_ss/sigma2[i]);
            double lambda_x_sum = 0.0;
            for(int j = 0; j<n; j++) {
                lambda_x_sum += x_arma(i,j)*lambda[j];
                if (isnan(lambda[j])) {throw std::runtime_error("lambda[j] should not be NaN");}
                if (isnan(x_arma(i,j))) {throw std::runtime_error("x_arma(i,j) should not be NaN");}
            }
            double m_alpha_0 = v_alpha_0*(lambda_x_sum/sigma2[i]);
            
            double v_alpha_1 = 1/(1/omega_1 + lambda_ss/sigma2[i]);
            double m_alpha_1 = v_alpha_1*(lambda_x_sum/sigma2[i]);

            if(z[i]) {
                alpha[i] = R::rnorm(m_alpha_0, sqrt(v_alpha_0));
            } else {
                alpha[i] = R::rnorm(m_alpha_1, sqrt(v_alpha_1));
            }

            double frac_1 = R::dnorm(0, m_alpha_0, sqrt(v_alpha_0), true) -
                R::dnorm(0, 0, sqrt(omega), true);
            double frac_2 = R::dnorm(0, 0, sqrt(omega_1), true) - 
                R::dnorm(0, m_alpha_1, sqrt(v_alpha_1), true);

            p_star[i] = p/(p + exp(frac_1 + frac_2)*(1-p));

            dinv(i,i) = 1/sigma2[i];

            if (sigma2[i]<0) {throw std::runtime_error("sigma2 should not be negative");}
            if (isnan(lambda_x_sum)) {throw std::runtime_error("lambda_x_sum should not be NaN");}
            if (v_alpha_0<0) {throw std::runtime_error("v_alpha_0 should not be negative");}
            if (v_alpha_1<0) {throw std::runtime_error("v_alpha_1 should not be negative");}
            if (isnan(m_alpha_0)) {throw std::runtime_error("m_alpha_0 should not be NaN");}
            if (isnan(m_alpha_1)) {throw std::runtime_error("m_alpha_1 should not be NaN");}
            if (isnan(sigma2[i])) {throw std::runtime_error("sigma should not be NaN");}
            if (isnan(alpha[i])) {throw std::runtime_error("alpha should not be NaN");}
            if (isnan(p_star[i])) {throw std::runtime_error("p_star should not be NaN");}
            if (isnan(dinv(i,i))) {throw std::runtime_error("dinv should not be NaN");}

            debug[i] = dinv(i,i);
        }

        for(int j = 0; j<n; j++) {
            arma::mat aux_1 = alpha.t()*dinv*alpha + 1;
            arma::mat aux_2 = alpha.t()*dinv*x_arma.col(j);
            if (aux_1.n_rows != 1) {throw std::runtime_error("aux_1 should have 1 row");}
            if (aux_1.n_cols != 1) {throw std::runtime_error("aux_1 should have 1 col");}
            if (!alpha.is_finite()) {throw std::runtime_error("alpha should be finite");}
            if (!dinv.is_finite()) {throw std::runtime_error("dinv should be finite");}
            if (!aux_1.is_finite()) {throw std::runtime_error("aux_1 should be finite");}
            

            double v_lambda = 1/(arma::as_scalar(aux_1));
            double m_lambda = v_lambda*(arma::as_scalar(aux_2));

            if (v_lambda<0) {throw std::runtime_error("v_lambda cannot be negative");}
            if (isnan(v_lambda)) {throw std::runtime_error("v_lambda should not be NaN");}
            if (isnan(m_lambda)) {throw std::runtime_error("m_lambda should not be NaN");}
            if (isnan(lambda[j])) {throw std::runtime_error("lambda should not be NaN");}

            lambda[j] = R::rnorm(m_lambda, sqrt(v_lambda));
        }

        alpha_matrix.row(k) = alpha;
        lambda_matrix.row(k) = lambda;
        sigma_matrix.row(k) = sigma2;
        z_matrix.row(k) = z;
        p_matrix.row(k) = p_star;
        debug_matrix.row(k) = debug;

    }
    return wrap(debug_matrix);             // Return to R
}
