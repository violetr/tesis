# SIMU 3.7.2
# normales multivariadas, estimacion de la matriz de precision
# grafico error vs n

#forma de grilla, distintos p, n aumentando, nrep repeticiones

#normales multivariadas

packages = c('MASS')

for (package in packages) {
  if (! require(package, character.only = T, quietly = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}


source('./metodos/seleccion_modelo_normal.R')
source('./auxiliares.R')

#################################################################

p = 120
nrep = 50

Theta = diag(p)

for (i in 1:p) {
    for (j in 1:p) {
    if (abs(j - i) == 1) {
      Theta[i, j] = 0.2
    }
  }
}

Sigma = solve(Theta)
mu = rep(0, p)

grilla <- floor(seq(50, 6000, length.out = 20))

error2 = rep(0, length(grilla))

for (iter in 1:nrep) {
  muestra = mvrnorm(N, mu, Sigma)
  
  for (k in 1:length(grilla)) {
    
    N = grilla[k]
    lambda = 2 * sqrt(log(p) /N)
    s = var(muestra)
    thestimada = glasso(s, rho = lambda, penalize.diagonal = FALSE)$wi
    error2[k] = error2[k] + base::norm(Theta - thestimada, "2")
  }
}

p = c(64, 120, 350)

pdf('./resultados/3_7_2_simu.pdf')
plot(grilla, error1 / (nrep), type = 'b', col = "blue", ylab = "", xlab = "")
legend("topright", 
       c("p=64","p=120","p=350"),
       lwd = 1, lty = c(1, 1), 
       pch = 1, 
       cex = 1.2, pt.cex = 1.5,
       col=c('blue', 'green', 'red'))
points(grilla, error2 / (nrep), type = 'b', col = "green")
points(grilla, error3 / (nrep), type = 'b', col = "red")
dev.off()
