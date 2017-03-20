# SIMU 3.2.1.3

######################################################################
#####analisis sesgo varianza lasso distintos valores de lambda########

set.seed(2017)

datotest=rnorm(1000)
ytest=datotest[1]*2+datotest[400]+datotest[3]*0.5+rnorm(1)
conjlambda=seq(0,1.7,0.05)
nlambda=length(conjlambda)

matrizestimacionesy=matrix(rep(0,nlambda*100),ncol=nlambda)

for(j in 1:100){
  
  datos=matrix(nrow=n,ncol=p)
  for( i in 1:p){
    datos[,i]=rnorm(n)
  }
  
  y=datos[,1]*2+datos[,400]+datos[,3]*0.5+rnorm(50,0,1)
  
  lassosv<-glmnet(x=datos,lambda=conjlambda,alpha=1,y=y)
  
  #matriz con coef filas y lambda como columnas
  
  conjlambda=lassosv$lambda
  
  ytestestimados=predict(lassosv,matrix(datotest,nrow=1),type="response")
  
  matrizestimacionesy[j,]=ytestestimados[1,]  
  
}

#calculo
VAR=apply(matrizestimacionesy,2,var)
MSE=apply((ytest-matrizestimacionesy)^2,2,mean)
SESGO=ytest-apply(matrizestimacionesy,2,mean)
#en un mismo grafico para cada lambda grafico estas tres funciones con lineas.  


#grafico
#paleta de colores
pal3 <- brewer.pal(10, "Set3")

plot(conjlambda,SESGO^2 ,  xlab = "",ylab="",type="l",col=pal3[3], xlim=c(0,1.8),ylim=c(0,9),lwd=2)
lines(conjlambda,MSE,col=pal3[1],lwd=2)
lines(conjlambda,VAR,col=pal3[6],lwd=2)
elegido=conjlambda[which.min(MSE)]
points(elegido,MSE[which.min(MSE)],pch=4,col=pal3[10],cex=2)


