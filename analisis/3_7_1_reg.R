# SIMU 3.7.1

# simulacion regresion altas dimensiones para aplicar lasso

# en la primera parte hago la comparacion entre lasso
# adaptive lasso,thresholded lasso y ridge para una simulacion 
# con variables normales. calculo el error cuadratico medio

# en la segunda parte hago un analisis del sesgo y la varianza para 
# lasso con distintos valores de lambda

#semilla
set.seed(2017)

source('./metodos/regresion_lineal.R')
source('./auxiliares.R')

###################################################################
#auxiliares


genero.muestra<-function(n,ntest,p){
  
  # genero matriz de datos
  datos=matrix(nrow=n,ncol=p)
  for( i in 1:p){
    datos[,i]=rnorm(n)
  }
  
  #genero datos para testear eficiencia
  datostest=matrix(nrow=ntest,ncol=p)
  for( i in 1:p){
    datostest[,i]=rnorm(ntest)
  }
  
  #variable respuesta con coef 1:2;3:1;400:0.5 y ruido normal con sigma 0.5
  y=datos[,1]*2+datos[,400]+datos[,3]*0.7+rnorm(n,0,0.5)
  ytest=datostest[,1]*2+datostest[,400]+datostest[,3]*0.7+rnorm(ntest,0,0.5)
  
  return(list(datostrain=datos,datostest=datostest,y=y,ytest=ytest))  
}

#####################################################################
############################# datos #################################

#numero de muestra y numero de variables
n=50
ntest=20
p=1000

# genero matriz de datos
datos=matrix(nrow=n,ncol=p)
for( i in 1:p){
  datos[,i]=rnorm(n)
}

#genero datos para testear eficiencia
datostest=matrix(nrow=ntest,ncol=p)
for( i in 1:p){
  datostest[,i]=rnorm(ntest)
}

#variable respuesta con coef 1:2;3:0.7;400:1 y ruido normal con sigma 0.5
y=datos[,1]*2+datos[,400]+datos[,3]*0.7+rnorm(n,0,0.5)
ytest=datostest[,1]*2+datostest[,400]+datostest[,3]*0.7+rnorm(ntest,0,0.5)

#####################################################################
############################ estimacion #############################

# voy a estimar los p coeff betas por lasso con 
# 10 folds CV para estimar lambda, uso paquete glmnet
# es un paquete que tiene un estimador elastic net 
# para regresiones lineales generalizadas
# alpha es el parametro de la combinacion
# convexa entre convexaridge y lasso, alpha=1 es lasso puro

##%%--=============== LASSO ===============--%%##

#ajusto modelo
modelolasso<-cv.glmnet(x=datos,y=y,alpha=1,nfolds = 10,intercept=FALSE)
#obtengo coeficientes para el lambda con menos error
coeflasso=as.matrix(coef(modelolasso,s="lambda.min"))[,1]
#saco el intercept
coeflasso=coeflasso[2:1001]

#coeficientes de interes
coeflasso[c(1,3,400)]

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2pathlasso.pdf")
#plot camino de coeficientes para distintos lambdas
plot(modelolasso$glmnet.fit,ylab="coeficientes beta", xvar = "lambda")
dev.off()
#calculo MSE para datos test
lassopred=predict(modelolasso, newx = datostest, s = "lambda.min",type = "response")
MSElasso=mean((lassopred-ytest)^2)


#cantidad de variables con coef distinto de cero
nnoceroslasso=p-sum(coeflasso==0)

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2lassocoef.pdf")
#grafico los coeficientes calculados y los verdaderos distintos de cero
plot(1:p,coeflasso,xlab="",ylab="coeficientes beta",ylim=c(-1,2.2))
points(c(1,3,400),c(2,0.7,1),col='red',pch=17)
dev.off()

##%%--=============== RIDGE ===============--%%##
#hago todo lo mismo pero con ridge

modeloridge<-cv.glmnet(x=datos,y=y,alpha=0,nfolds = 10,intercept=FALSE)
coefridge=as.matrix(coef(modeloridge,s="lambda.min"))[,1]
coefridge=coefridge[2:1001]

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2pathridge.pdf")
plot(modeloridge$glmnet.fit, xvar = "lambda",ylab="coeficientes beta")
dev.off()

ridgepred=predict(modeloridge, newx = datostest, s = "lambda.min")
MSEridge=mean((ridgepred-ytest)^2)

nnocerosridge=p-sum(coefridge==0)

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2ridgecoef.pdf")
plot(1:p,coefridge,xlab="",ylab="coeficientes beta",ylim=c(-1,2.2))
points(c(1,3,400),c(2,0.7,1),col='red',pch=17)
dev.off()

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2ridgecoef2.pdf")
plot(1:p,coefridge,xlab="",ylab="coeficientes beta")
points(c(1,3,400),c(2,0.7,1),col='red',pch=17)
dev.off()

##%%--=========== ADAPTIVE LASSO (dos pasos) ===========--%%##

#Genero pesos para que los coeficientes chicos de lasso
#tengan penalidad muymuy alta, ajusto lasso con esos pesos

pesos=as.integer(!coeflasso==0)*1/abs(coeflasso)
pesos[is.na(pesos)]=10000

modeloadaplasso<-cv.glmnet(x=datos,y=y,alpha=1,nfolds = 10,penalty.factor = pesos,intercept=FALSE)

#hago el mismo analisis que en los dos metodos anteriores

coefadap=as.matrix(coef(modeloadaplasso,s="lambda.min"))[,1]
coefadap=coefadap[2:1001]

plot(modeloadaplasso$glmnet.fit,xvar="lambda")

adaplassopred=predict(modeloadaplasso, newx = datostest, s = "lambda.min")


MSEadaplasso=mean((adaplassopred-ytest)^2)

nnocerosadap=p-sum(coefadap==0)

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2adaplassocoef.pdf")
plot(1:p,coefadap,xlab="",ylab="coeficientes beta",ylim=c(-1,2.2))
points(c(1,3,400),c(2,0.7,1),col='red',pch=17)
dev.off()

##%%--============= THRESHOLDED LASSO =============--%%##
# tengo que elegir el valor limite tau

###############

#hago el mismo analisis que para los demas metodos

grupos=kfoldgrupos(n,5)
modelothres=cv.tau.threslasso(coeflasso,modelolasso$lambda.min,grupos,datos,y,100)

coefthres=modelothres$coef
tau=modelothres$tau

threslassopred=predecir(datostest,coefthres)

MSEthreslasso=mean((t(threslassopred)-ytest)^2)
p=1000
nnocerosthres=p-sum(coefthres==0)

pdf("C:/Users/violeta/Dropbox/tesisvio/cap2threslassocoef.pdf")
plot(1:p,coefthres,xlab="",ylab="coeficientes beta",ylim=c(-1,2.2))
points(c(1,3,400),c(2,0.7,1),col='red',pch=17)
dev.off()


######################################################################
# SIMU nrep
# simulacion para medir MSE

n=50
ntest=20
p=1000
nrep=20


set.seed(2017)

nrep=100

n=50
p=1000
ntest=20
MSE=rep(0,4)
nnoceros=rep(0,4)
incluyetruevar=rep(0,4)

for(iter in 1:nrep){
  
  muestra=genero.muestra(n,ntest,p)
  datos=muestra$datostrain
  y=muestra$y
  datostest=muestra$datostest
  ytest=muestra$ytest
  
  modelolasso<-cv.glmnet(x=datos,y=y,alpha=1,nfolds = 5,intercept=FALSE)
  #obtengo coeficientes para el lambda con menos error
  coeflasso=coef(modelolasso,s="lambda.min")
  #saco el intercept
  length(coeflasso)
  hasta=p+1
  coeflasso=coeflasso[2:hasta]

  lassopred=predict(modelolasso, newx = datostest, s = "lambda.min",type = "response")
  lassopred==predecir(datostest,coeflasso)
  
  MSE[1]=MSE[1]+mean((lassopred-ytest)^2)
  
  nnoceros[1]=nnoceros[1]+p-sum(coeflasso==0)
  
  incluyetruevar[1]=incluyetruevar[1]+(coeflasso[1]!=0 & coeflasso[3]!=0 & coeflasso[400]!=0)
  
  ##%%--=============== RIDGE ===============--%%##
  #hago todo lo mismo pero con ridge
  
  modeloridge<-cv.glmnet(x=datos,y=y,alpha=0,nfolds = 5,intercept=FALSE)
  coefridge=as.matrix(coef(modeloridge,s="lambda.min"))[,1]
  coefridge=coefridge[2:hasta]
  
  ridgepred=predict(modeloridge, newx = datostest, s = "lambda.min")
  ridgepred==predecir(datostest,coefridge)
  MSE[2]=MSE[2]+mean((ridgepred-ytest)^2)
  
  nnoceros[2]=nnoceros[2]+p-sum(coefridge==0)
  
  incluyetruevar[2]=incluyetruevar[2]+(coefridge[1]!=0 & coefridge[3]!=0 & coefridge[400]!=0)
  
  ##%%--=========== ADAPTIVE LASSO (dos pasos) ===========--%%##
  
  pesos=as.integer(!coeflasso==0)*1/abs(coeflasso)
  pesos[is.na(pesos)]=10000
  
  modeloadaplasso<-cv.glmnet(x=datos,y=y,alpha=1,nfolds = 10,penalty.factor = pesos,intercept=FALSE)
  
  coefadap=coef(modeloadaplasso,s="lambda.min")
  coefadap=coefadap[2:hasta]

  adaplassopred=predict(modeloadaplasso, newx = datostest, s = "lambda.min")
  adaplassopred==predecir(datostest,coefadap)
  MSE[3]=MSE[3]+mean((adaplassopred-ytest)^2)
  
  nnoceros[3]=nnoceros[3]+p-sum(coefadap==0)
  
  incluyetruevar[3]=incluyetruevar[3]+(coefadap[1]!=0 & coefadap[3]!=0 & coefadap[400]!=0)
  
  #thres
  
  grupos=kfoldgrupos(n,5)
  modelothres=cv.tau.threslasso(coeflasso,modelolasso$lambda.min,grupos,datos,y,100)
  
  coefthres=modelothres$coef
  taucv[iter]=modelothres$tau
  
  threslassopred=predecir(datostest,coefthres)
  
  MSE[4]=MSE[4]+mean((t(threslassopred)-ytest)^2)
  
  nnoceros[4]=nnoceros[4]+p-sum(coefthres==0)

  incluyetruevar[4]=incluyetruevar[4]+(coefthres[1]!=0 & coefthres[3]!=0 & coefthres[400]!=0)

}

MSE=MSE/nrep
nnoceros= nnoceros/nrep
incluyetruevar=incluyetruevar/nrep


