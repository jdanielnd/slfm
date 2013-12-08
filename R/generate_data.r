sample_data <- function(m = 40, n = 251) {
	m = 40
	n = 251

	lambda = array(0,c(1,n))
	for(j in 1:n){lambda[j] = rnorm(1,0,1)}
	alpha = array(0,c(m,1))
	for(i in 1:m){alpha[i] = rnorm(1,0,1)}
	alpha.zero = sample(1:m,35)
	alpha[alpha.zero] = 0
	sigma2 = array(0,c(m,1))
	for(i in 1:m){sigma2[i] = 1/rgamma(1, 2.1,rate=0.55)}
	epsilon = array(0,c(m,n))
	for(i in 1:m){epsilon[i,] = rnorm(n,0,sqrt(sigma2[i]))}

	# MATRIZ X
	x = alpha%*%lambda+epsilon
	x
}