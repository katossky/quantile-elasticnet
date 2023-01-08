clear()

########################### Simulation des données #############################

## choix des paramètres
n = 500
tau0 <- 2
beta0 <- c(3, 3, 0, 3, 3)
l <- length(beta0)

## quantile 
theta <- 0.5
xi1 <- (1 - 2*theta) / (theta * (1-theta))
xi2 <- sqrt(2 / (theta * (1-theta)))

## simulation des variables v, z, u, x, et y
v <- rexp(n)
v_bar <- v / tau0
z <- rnorm(n)
u <- xi1 * v_bar + xi2 * sqrt(v_bar / tau0) * z
x <- array(rnorm(n*l), c(l, n))
y <- t(x) %*% beta0 + u

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
eta2 <- rgamma(1, sgape=c2, rate=d2)
beta <- rep(0, l)



#x est un vecteur contenant une série de valeurs
cv<-function(x){
  moy<-mean(x) # moyenne de x
  s<-sd(x)# ecart type de x
  rslt<-s/moy # calcul du CV
  rslt #la fonction retourne le résultat
}