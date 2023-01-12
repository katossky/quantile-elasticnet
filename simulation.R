rm(list=ls())

library('statmod')
library('GIGrvg')
library('mvtnorm')
library('hqreg')

########################### Génération des données #############################

# fonctions Xi
Xi1 <- function(theta){
  return ( (1 - 2*theta) / (theta * (1-theta)) )
}

Xi2 <- function(theta){
  return (sqrt(2 / (theta * (1-theta))))
}

data_generation <- function(n, theta, beta0=c(3, 3, 0, 3, 3), tau0=2, g_err=TRUE){
  # renvoie les un tableau dont la première colonne y, les colonnes suivantes x
  # n nombre d échantillons à générer (=nombre de lignes), 
  # theta quantile souhaité
  
  xi1 <- Xi1(theta)
  xi2 <- Xi2(theta)
  p <- length(beta0)
  
  # simulation des variables explicatives x
  S <- matrix(0, nrow=p, ncol=p)
  for(k in 1:p){
    for(l in 1:p){
      S[k,l] <- 0.5**abs(k-l)
    }
  }
  
  x <- matrix(0, nrow=n, ncol=p)
  for(i in 1:n){
    x[i,] <- rmvnorm(1, mean=rep(0, p), sigma=S)
  }
  
  # simulation des erreurs 
  
  if(g_err==FALSE){
  ## erreurs du modèle quantile elastic net :
  v <- rexp(n)
  v_bar <- v / tau0
  z <- rnorm(n)
  u <- xi1 * v_bar + xi2 * sqrt(v_bar / tau0) * z
  }
  
  if(g_err){
  ## erreurs gaussiennes :
  delta <- qnorm(theta, mean=0, sd=3)
  u <- rnorm(n, mean=-delta, sd=3)
  }
  
  # calcul de y :
  y <- x %*% beta0 + u
  
  list(x=x, y=y)
}


########################### Fonction Update ###############################

a <- 0.1
b <- 0.1
c1 <- 0.1
c2 <- 0.1
d1 <- 0.1
d2 <- 0.1

v_bar_update <- function(tau, beta, n, x, y, xi2, xi1){
  chi_v <- tau * (y - x%*%beta)**2 / (xi2 ** 2)
  psi_v <-(tau * xi1**2)/xi2**2 + 2*tau 
  
  v_bar <- rep(0,n)
  for(i in 1:n){
  v_bar[i] <- rgig(1, lambda=1/2, chi=chi_v[i], psi=psi_v)
  }
  return(v_bar)
}

y_bar_update <- function(v_bar, beta, n, p, xi1, x, y){
  y_bar <- matrix(nrow=n, ncol=p)
  for(k in 1:p){
    y_bar[,k] <- y - xi1*v_bar - rowSums( t( t(x)*beta)[,-k]  )
  }
  
  return(y_bar)
}

t_update <- function(beta, eta1_bar, eta2, p){
  chi_t <- 2 * eta2 * (beta**2)
  psi_t <- 2 * eta1_bar
  
  t <- rep(0, p)
  for(k in 1:p){
    t[k] <- 1 + rgig(1, lambda=1/2, chi=chi_t[k], psi=psi_t)
  }
  return(t)
}

beta_update <- function(tau, v_bar, t, y_bar, eta2, xi2, x, y, p){
  
  sigma_inv2 <- tau * xi2**(-2) * colSums((x**2) * v_bar**(-1)) + 2 * eta2 * t*(t - 1)**(-1)
  sigma_bar <- 1 / sqrt(sigma_inv2)
  mu_bar <- (sigma_bar**2) * tau * (xi2**(-2)) * colSums( y_bar * x * (v_bar**(-1)) )
  
  beta <- rnorm(p, mean=mu_bar, sd=sigma_bar)
  
  return(beta)
}

eta1_bar_update <- function(t, eta1_bar, p){
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

eta2_update <- function(beta, t, p){
  eta2 <- rgamma( 1, shape=c2 + p/2, rate=d2 + sum( (beta**2)*t/(t-1) ) )
  return(eta2)
}

tau_update <- function(v_bar, beta, x, y, xi1, xi2, n){
  r <- b + sum( v_bar + (y - x%*%beta - xi1*v_bar)**2 / (2*v_bar*xi2**2) )
  tau <- rgamma( 1, shape=a + 3*n/2, rate= r)
  return(tau)
}

Gibbs_update <- function(tau, eta1_bar, eta2, beta, x, y, theta){
  # mets à jours les paramètres tau, eta1_bar, eta2 et beta après une étape de Gibbs
  
    n <- length(y)
    p <- length(beta)
    xi1 <- Xi1(theta)
    xi2 <- Xi2(theta)
  
    v_bar <- v_bar_update(tau, beta, n, x, y, xi2, xi1)
    y_bar <- y_bar_update(v_bar, beta, n, p, xi1, x, y)
    t <- t_update(beta, eta1_bar, eta2, p)
    new_beta <- beta_update(tau, v_bar, t, y_bar, eta2, xi2, x, y, p)
    new_eta1_bar <- eta1_bar_update(t, eta1_bar, p)
    new_eta2 <- eta2_update(new_beta, t, p)
    new_tau <- tau_update(v_bar, new_beta, x, y, xi1, xi2, n)
    
    return( c(new_tau, new_eta1_bar, new_eta2, new_beta) )
}

########################### Simulation ###############################

# pour un échantillon iid prendre 
#           burnin = 10 
#           autocorr = 10

Generation_Gibbs <- function(r_gibbs, x, y, theta, burnin=0, autocorr=0){
  
  # génère un tableau avec r_gibbs échantillonages de Gibbs des variables :
  # tau : première colonne
  # eta1_bar : deuxième colonne
  # eta2 : troisième colonne
  # beta : colonnes restantes

  # initialisation des paramètres 
  tau <- rgamma(1, shape=a, rate=b)
  eta1 <- rgamma(1, shape=c1, rate=d1)
  eta2 <- rgamma(1, shape=c2, rate=d2)
  eta1_bar <- (eta1**2) / (4*eta2)
  beta <- rep(1, length(x[1,]))

  # tables vides
  tau_gibbs <- matrix(nrow=r_gibbs, ncol=1)
  eta1_bar_gibbs <- matrix(nrow=r_gibbs, ncol=1)
  eta2_gibbs <- matrix(nrow=r_gibbs, ncol=1)
  beta_gibbs <- matrix(nrow=r_gibbs, ncol=length(x[1,]))

  # burn in :
  for(r in 1:burnin){
    res <- Gibbs_update(tau, eta1_bar, eta2, beta, x, y, theta)
    
    tau <- res[1]
    eta1_bar <- res[2]
    eta2 <- res[3]
    beta <- res[-c(1:3)]
  }
  
  # gibbs sampling :
  for(r in 1:r_gibbs){
      for(j in 1:(1+autocorr)){
        res <- Gibbs_update(tau, eta1_bar, eta2, beta, x, y, theta)
    
        tau <- res[1]
        eta1_bar <- res[2]
        eta2 <- res[3]
        beta <- res[-c(1:3)]
      }
    
      tau_gibbs[r] <- tau
      eta1_bar_gibbs[r] <- eta1_bar
      eta2_gibbs[r] <- eta2
      beta_gibbs[r,] <- beta
  }
  
  #return(cbind(tau_gibbs, eta1_bar_gibbs, eta2_gibbs, beta_gibbs))
  list(tau = tau_gibbs, eta1_bar = eta1_bar_gibbs, eta2 = eta2_gibbs, beta = beta_gibbs)

}

bayesian_qreg <- function(x, y, theta, n_vals=20){
  # genere un estimateur bayesien de beta en moyennant sur n_vals valeurs generees par gibbs sampling
  # output : 
  # premiere ligne estimateur de beta
  # deuxième ligne : standard deviation
  
  gen <- Generation_Gibbs(n_vals, x, y, theta, burnin=10, autocorr=10)
  mean <- colMeans(gen$beta)
  std <- sqrt( colMeans(t(t(gen$beta) - mean)**2) )
  
  l2 <- mean(gen$eta2 / gen$tau)
  l1 <- mean(sqrt(4 * gen$eta2 * gen$eta1_bar) / gen$tau)
  
  list(beta_mean=mean, beta_std=std, tau=mean(gen$tau), lambda1=l1, lambda2=l2)
}

################## paramétrisation burn-in & corrélations ######################

theta = 0.5
data <- data_generation(100, theta)
y <- data$y
x <- data$x

r_gibbs = 200

res <- Generation_Gibbs(r_gibbs, x, y, theta)
tau_gibbs <- res$tau
eta1_bar_gibbs <- res$eta1_bar
eta2_gibbs <- res$eta2
beta_gibbs <- res$beta

# auto-corrélations / burn-in ?
plot(1:30,beta_gibbs[1:30,1]) 
plot(1:30,beta_gibbs[1:30,3])
plot(1:30,tau_gibbs[1:30]) 
plot(1:200,eta1_bar_gibbs[1:200]) 
plot(1:200,eta2_gibbs[1:200]) 
  # --> burn-in : enlever 10 premières valeurs 

acf(beta_gibbs[,1])
acf(beta_gibbs[,3])
acf(tau_gibbs)
acf(eta1_bar_gibbs)
acf(eta2_gibbs)
  # --> correlations : prendre une valeur toutes les 5 ou 10 valeurs 

rm(beta_gibbs, r_gibbs, data, eta1_bar_gibbs, eta2_gibbs, res, tau_gibbs, x, y, theta)

################################ Analyse #################################

theta = 0.5
data <- data_generation(100, theta)
y <- data$y
x <- data$x

r_gibbs = 100
burnin = 10
autocorr = 10

res <- Generation_Gibbs(r_gibbs, x, y, theta, burnin, autocorr)
tau_gibbs <- res$tau
eta1_bar_gibbs <- res$eta1_bar
eta2_gibbs <- res$eta2
beta_gibbs <- res$beta

hist(tau_gibbs)
hist(eta1_bar_gibbs)
hist(eta2_gibbs)
hist(beta_gibbs[,1])
hist(beta_gibbs[,3])

# estimateurs :
sprintf('variable tau : moyenne %f, écart-type %f', mean(tau_gibbs), sd(tau_gibbs))
sprintf('variable eta1 : moyenne %f, écart-type %f', mean(eta1_bar_gibbs), sd(eta1_bar_gibbs))
sprintf('variable eta2 : moyenne %f, écart-type %f', mean(eta2_gibbs), sd(eta2_gibbs))

sprintf('variable beta1 : moyenne %f, écart-type %f', mean(beta_gibbs[,1]), sd(beta_gibbs[,1]))
sprintf('variable beta2 : moyenne %f, écart-type %f', mean(beta_gibbs[,2]), sd(beta_gibbs[,2]))
sprintf('variable beta3 : moyenne %f, écart-type %f', mean(beta_gibbs[,3]), sd(beta_gibbs[,3]))

rm(beta_gibbs, r_gibbs, data, eta1_bar_gibbs, eta2_gibbs, res, tau_gibbs, x, y, theta, autocorr, burnin)

######################### Estimateur avec hqreg ##########################

classical_qreg <- function(x, y, theta){
  alpha <- c(0:16)*0.05 + rep(0.1, 17)
  
  best_error <- 10**50
  best_alpha <- 0
  best_lambda <- 0
  best_beta <- rep(0, length(x[1,]))
  
  for(i in 1:length(alpha)){
    reg <- cv.hqreg(x, y, alpha=alpha[i], tau=theta, lambda.min=0.01, method='quantile', FUN='hqreg', type.measure='mse', seed=42)
    s <- reg$lambda == reg$lambda.min
    error <- sum(s*reg$cve)
    
    print(reg$lambda.min)
    
      if(error < best_error){
        best_error <- error
        best_alpha <- alpha[i]
        best_lambda <- reg$lambda.min
        best_beta <- coef(reg, lambda=best_lambda)
      }
  }
    list(alpha=best_alpha, lambda=best_lambda, beta=best_beta[-1])
}

######################### Fonction perte ##################################

quant_loss <- function(z, theta){
  if(z >= 0){
    l = theta*z
  }
  else{
    l = - (1 - theta)*z
  }
  return(l)
}

Elastic_net_loss <- function(x, y, theta, beta, l1=0, l2=0, lambda=0, alpha=0, type='bayesian'){
  n = length(y)
  
  loss = 0
  for (i in 1:n){
    loss = loss + quant_loss(y[i] - sum(x[i,]*beta), theta) / n
  }
  
  if(type=='bayesian'){
    loss <- loss + l1*sum(beta) + l2*sum(beta*beta)
  }
  
  else{
    loss <- loss + lambda * alpha * sum(beta) + lambda * (1-alpha) * sum(beta*beta) / 2
  }

  return(loss)  
}

############################ Comparaisons ################################

n_samples_train = 50
n_samples_test = 100
theta = 0.3
beta0 = c(rep(3,5), rep(0, 3), rep(3, 5))
tau0 = 2
g_err = FALSE


training_data <- data_generation(n_samples_train, theta, beta0,  tau0, g_err)
x_train <- training_data$x
y_train <- training_data$y

testing_data <- data_generation(n_samples_test, theta, beta0, tau0, g_err)
x_test <- testing_data$x
y_test <- testing_data$y

# classical regression
t1 <- Sys.time()
qreg <- classical_qreg(x_train, y_train, theta)
t2 <- Sys.time()
print(t2 - t1)
c_loss <- Elastic_net_loss(x_test, y_test, theta, qreg$beta, lambda=qreg$lambda, alpha=qreg$alpha, type='classical')
print(c_loss)
   # note for g_err = FALSE / theta = 0.3
   # elapsed time : 39s
   # loss : 0.7

# bayesian regression
n_values = c(100, 200, 300, 400, 500)
b_loss = rep(0,5)
for(l in 1:5){
    print('nb d iterations de gibbs :')
    print(n_values[l])
    t1 <- Sys.time()
    bayes_qreg <- bayesian_qreg(x_train, y_train, theta, n_vals=n_values[l])
    t2 <- Sys.time()
    print(t2 - t1)
    b_loss[l] <- Elastic_net_loss(x_test, y_test, theta, bayes_qreg$beta_mean, l1=bayes_qreg$lambda1, l2=bayes_qreg$lambda2)
    print(b_loss[l])
}
  # note for g_err = FALSE / theta = 0.3
  # elapsed time : ~ 1.5 s pour 100 gibbs steps
  # loss : toujours autour de 9

sum((bayes_qreg$beta_mean - beta0)**2) 
sum((qreg$beta - beta0)**2)
# --> la méthode bayésienne est plus proche !
# revoir calcul des lambdas
