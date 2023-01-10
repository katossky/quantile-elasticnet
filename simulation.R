rm()

########################### Simulation des données #############################

## choix des paramètres
n = 500
tau <- 2
beta0 <- c(3, 3, 0, 3, 3)
p <- length(beta0)

## quantile 
theta <- 0.5
xi1 <- (1 - 2*theta) / (theta * (1-theta))
xi2 <- sqrt(2 / (theta * (1-theta)))

## simulation des variables v, z, u, x, et y
v <- rexp(n)
v_bar <- v / tau
z <- rnorm(n)
u <- xi1 * v_bar + xi2 * sqrt(v_bar / tau) * z

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

r_gibbs <- 1000

Gibbs_update <- function(x, y, tau, xi1, xi2, eta1_bar, eta2, beta){
    # v_bar 
    lambda <- tau * (y - x%*%beta)**2 / (xi2 ** 2)
    mu2 <- lambda / ( (tau * xi1**2)/xi2**2 + 2*tau ) 
    mu <- sqrt(mu2)
    v_bar <- rinvgauss(1, mean = mu, shape = lambda)
    
    # y_bar 
    y_bar <- matrix(nrow=n, ncol=p)
    for(k in 1:p){
        y_bar[,k] <- y - xi1*v_bar - sum( t( (t(x)*beta)[-k] ) )
    }
    
    # t
    t <- rep(2, p) ### A completer
    
    # sigma
    #sigma_inv2[k] <- tau * xi2**(-2) * sum(x[,k] * v_bar**(-1)) + 2 * eta2 * t[k]*(t[k] - 1)**(-1)
    sigma_inv2 <- tau * xi2**(-2) * colSums(x * v_bar**(-1)) + 2 * eta2 * t*(t - 1)**(-1)
  
    eta1_bar; eta2; beta
}
