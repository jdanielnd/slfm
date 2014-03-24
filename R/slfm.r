slfm <- function(
  x, ite, a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1,
  omega = 10, omega_1 = 0.01, burnin = 500) {
  
  res <- bfmu_c(x, ite, a, b, gamma_a, gamma_b, omega, omega_1)
  coda::as.mcmc(res[burnin:ite,])
}