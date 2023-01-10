rm()

library('statmod')

########################### Simulation des données #############################

## choix des paramètres
n = 500
tau0 <- 2
beta0 <- c(3, 3, 0, 3, 3)
p <- length(beta0)

## quantile 
theta <- 0.5
xi1 <- (1 - 2*theta) / (theta * (1-theta))
xi2 <- sqrt(2 / (theta * (1-theta)))

## simulation des variables v, z, u, x, et y
v <- rexp(n)
v_bar <- v / tau0
z <- rnorm(n)
u <- xi1 * v_bar + xi2 * sqrt(v_bar / tau0) * z

x <- matrix(rnorm(n*p), nrow=n, ncol=p)
y <- x %*% beta0 + u

########################### Lois Conditionnelles ###############################

a <- 1
b <- 1
c1 <- 1
c2 <- 1
d1 <- 1
d2 <- 1

# initialisation des paramètres 
tau <- rgamma(1, shape=a, rate=b)
eta1 <- rgamma(1, shape=c1, rate=d1)
eta2 <- rgamma(1, shape=c2, rate=d2)
eta1_bar <- (eta1**2) / (4*eta2)
beta <- rep(0, p)

Gibbs_update <- function(tau, eta1_bar, eta2, beta){
    # v_bar 
    lambda <- tau * (y - x%*%beta)**2 / (xi2 ** 2)
    mu2 <- lambda / ( (tau * xi1**2)/xi2**2 + 2*tau ) 
    mu <- sqrt(mu2)
    
    v_bar <- rep(0,n)
    for(i in 1:n){
      v_bar[i] <- rinvgauss(1, mean = mu[i], shape = lambda[i])
    }
    
    # y_bar
    y_bar <- matrix(nrow=n, ncol=p)
    for(k in 1:p){
        y_bar[,k] <- y - xi1*v_bar - sum( t( (t(x)*beta)[-k] ) )
    }
    
    # t
    t <- rep(2, p) ### A completer
    
    # Beta
    sigma_inv2 <- tau * xi2**(-2) * colSums((x**2) * v_bar**(-1)) + 2 * eta2 * t*(t - 1)**(-1)
    sigma_bar <- 1 / sqrt(sigma_inv2)
    # !! production de NAN
    mu_bar <- (sigma_bar**2) * tau * (xi2**(-1)) * colSums( y_bar * x * (v_bar**(-1)) )
    
    new_beta <- rnorm(p, mean=mu_bar, sd=sigma_bar)
    
    # eta1_bar
    new_eta1_bar <- rgamma(1, shape=c1 + p/2, rate=d1 + sum(t))
    
    # eta2
    new_eta2 <- rgamma( 1, shape=c2 + p/2, rate=d2 + sum( (new_beta**2)*t/(t-1) ) )
    
    # tau
    new_tau <- rgamma( 1, shape=a + 3*n/2, rate= b + sum( v_bar + (y - x%*%new_beta - xi1*v_bar)**2 / (2*(xi2**2)*v_bar) ))
    
    return( c(new_tau, new_eta1_bar, new_eta2, new_beta) )
}

res <- Gibbs_update(tau, eta1_bar, eta2, beta)
tau <- res[1]
eta1_bar <- res[2]
eta2 <- res[3]
beta <- res[-c(1:3)]

r_gibbs <- 1000


