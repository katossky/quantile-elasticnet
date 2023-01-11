rm()

library('statmod')
library('GIGrvg')

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

########################### Fonction Update ###############################

a <- 1
b <- 1
c1 <- 1
c2 <- 1
d1 <- 1
d2 <- 1

v_bar_update <- function(tau, beta){
  chi_v <- tau * (y - x%*%beta)**2 / (xi2 ** 2)
  psi_v <-(tau * xi1**2)/xi2**2 + 2*tau 
  
  v_bar <- rep(0,n)
  for(i in 1:n){
  v_bar[i] <- rgig(1, lambda=1/2, chi=chi_v[i], psi=psi_v)
  }
  return(v_bar)
}

y_bar_update <- function(v_bar, beta){
  y_bar <- matrix(nrow=n, ncol=p)
  for(k in 1:p){
    y_bar[,k] <- y - xi1*v_bar - rowSums( t( t(x)*beta)[,-k]  )
  }
  
  return(y_bar)
}

t_update <- function(beta, eta1_bar, eta2){
  chi_t <- 2 * eta2 * (beta**2)
  psi_t <- 2 * eta1_bar
  
  t <- rep(0, p)
  for(k in 1:p){
    t[k] <- 1 + rgig(1, lambda=1/2, chi=chi_t[k], psi=psi_t)
  }
  return(t)
}

beta_update <- function(tau, v_bar, t, y_bar, eta2){
  
  sigma_inv2 <- tau * xi2**(-2) * colSums((x**2) * v_bar**(-1)) + 2 * eta2 * t*(t - 1)**(-1)
  sigma_bar <- 1 / sqrt(sigma_inv2)
  mu_bar <- (sigma_bar**2) * tau * (xi2**(-2)) * colSums( y_bar * x * (v_bar**(-1)) )
  
  beta <- rnorm(p, mean=mu_bar, sd=sigma_bar)
  
  return(beta)
}

eta1_bar_update <- function(t, eta1_bar){
  # one step Metropolis Hastings
  # proposition pour new_eta1_bar :
  prop_eta1_bar <- rgamma(1, shape=c1 + p, rate=d1 + sum(t-1))
  
  # acceptation / rejet :
  alpha = dgamma(prop_eta1_bar, shape=c1 + p, rate=d1 + sum(t-1)) / dgamma(eta1_bar, shape=c1 + p, rate=d1 + sum(t-1))
  u = runif(1)
  
  new_eta1_bar = eta1_bar
  if(u<alpha){new_eta1_bar = prop_eta1_bar}
  return(new_eta1_bar)
}

eta2_update <- function(beta, t){
  eta2 <- rgamma( 1, shape=c2 + p/2, rate=d2 + sum( (beta**2)*t/(t-1) ) )
  return(eta2)
}

tau_update <- function(v_bar, beta){
  r <- b + sum( v_bar + (y - x%*%beta - xi1*v_bar)**2 / (2*v_bar*xi2**2) )
  tau <- rgamma( 1, shape=a + 3*n/2, rate= r)
  return(tau)
}


Gibbs_update <- function(tau, eta1_bar, eta2, beta){
    
    v_bar <- v_bar_update(tau, beta)
    y_bar <- y_bar_update(v_bar, beta)
    t <- t_update(beta, eta1_bar, eta2)
    new_beta <- beta_update(tau, v_bar, t, y_bar, eta2)
    new_eta1_bar <- eta1_bar_update(t, eta1_bar)
    new_eta2 <- eta2_update(new_beta, t)
    new_tau <- tau_update(v_bar, new_beta)
    
    return( c(new_tau, new_eta1_bar, new_eta2, new_beta) )
}

########################### Simulation ###############################

# initialisation des paramètres 
tau <- rgamma(1, shape=a, rate=b)
eta1 <- rgamma(1, shape=c1, rate=d1)
eta2 <- rgamma(1, shape=c2, rate=d2)
eta1_bar <- (eta1**2) / (4*eta2)
beta <- rep(1, p)

r_gibbs <- 100
tau_gibbs <- matrix(nrow=r_gibbs, ncol=1)
eta1_bar_gibbs <- matrix(nrow=r_gibbs, ncol=1)
eta2_gibbs <- matrix(nrow=r_gibbs, ncol=1)
beta_gibbs <- matrix(nrow=r_gibbs, ncol=p)

for(r in 1:r_gibbs){
    res <- Gibbs_update(tau, eta1_bar, eta2, beta)
    
    tau <- res[1]
    eta1_bar <- res[2]
    eta2 <- res[3]
    beta <- res[-c(1:3)]
    
    tau_gibbs[r] <- tau
    eta1_bar_gibbs[r] <- eta1_bar
    eta2_gibbs[r] <- eta2
    beta_gibbs[r,] <- beta
}

########################### Analyse ###############################

# auto-corrélations / burn-in ?
hist(beta_gibbs[,1])
plot(1:r_gibbs,beta_gibbs[1:r_gibbs,1])
