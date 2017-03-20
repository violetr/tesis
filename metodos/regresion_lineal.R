# METODOS reglineal

source('./auxiliares.R')


# regresion lineal:
#   adaptive lasso cv
#   thresholded lasso cv


#adaptive lasso

cv.adaptive.lasso<-function(datos,y,nfolds,termino.indep=TRUE){
  
  modelolasso<-cv.glmnet(x=datos,y=y,alpha=1,nfolds = 10,intercept=FALSE)
  #obtengo coeficientes para el lambda con menos error
  coeflasso=coef(modelolasso,s="lambda.min")
  #saco el intercept
  coeflasso=coeflasso[2:1001]
  
  pesos=as.integer(!coeflasso==0)*1/abs(coeflasso)
  pesos[is.na(pesos)]=10000
  
  modeloadaplasso<-cv.glmnet(x=datos,y=y,alpha=1,nfolds = nfolds,penalty.factor = pesos,intercept=termino.indep)
  
  #hago el mismo analisis que en los dos metodos anteriores
  
  coefadap=coef(modeloadaplasso,s="lambda.min")
  if(!intercept){
    coefadap=coefadap[2:1001]
  }
  
  return(coefadap)

}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# FUNCION
# threslasso(datosexplicativos,respuesta,tau,lambda):
# metodo thresholded lasso para regresion

# INPUT:
# datosexplicativos: matriz con datos explicativos dim=nx(p-1)
# respuesta: vector respuesta length=n
# tau: valor critico de corte
# lambda: valor de regularizacion

# OUTPUT:
# sthres: variables que fueron refiteadas
# coefthres: coeficientes de la regresion  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
threslasso<-function(datosexplicativos,respuesta,tau,lambda){
  
  modelolasso<-glmnet(x=datosexplicativos,y=respuesta,alpha=1,lambda = lambda)
  coeflasso=as.matrix(modelolasso$beta)[,1]
  
  coefthreslasso=coeflasso*as.integer(abs(coeflasso)>tau)
  seleccionada=(coefthreslasso !=0)
  
  sthres=obtenerindices(seleccionada)
  
  if(!(length(sthres) && sthres==0)){
    
    coefthresrefit=as.vector(t(matrix.inverse(t(datosexplicativos[,sthres])%*%t(t(datosexplicativos[,sthres])))%*%t(datosexplicativos[,sthres])%*%respuesta))
    for(k in 1:length(sthres)){
      coefthreslasso[sthres[k]]=coefthresrefit[k]  
    }
  }
  
  return(list(sthres=sthres,coefthres=coefthreslasso))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# FUNCION
# cv.tau.threslasso(datosexplicativos,respuesta,tau,lambda):
# metodo thresholded lasso para regresion 

# INPUT:
# beta.init: vector de coeficientes de la estimacion
# inicial con lasso
# lambda.best: lambda.cv de la estimacion inicial
# con lasso
# datosexplicativos: matriz con datos explicativos dim=nx(p-1)
# respuesta: vector respuesta length=n
# ntau: cantidad de 

# OUTPUT:
# sthres: variables que fueron refiteadas
# coefthres: coeficientes de la regresion  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


cv.threslasso<-function(coeflasso,datosexplicativos,respuesta,tau){
  
  coefthreslasso=coeflasso*as.integer(abs(coeflasso)>tau)
  seleccionada=(coefthreslasso !=0)
  
  sthres=obtenerindices(seleccionada)
  
  
  if(!(length(sthres) && sthres==0)){
    
    coefthresrefit=lm(respuesta ~ . - 1,as.data.frame(datosexplicativos[,sthres]))$coefficients
    
    for(k in 1:length(sthres)){
      coefthreslasso[sthres[k]]=coefthresrefit[k]  
    }
    
  }
  
  return(list(sthres=sthres,coefthres=coefthreslasso))
}


cv.tau.threslasso<-function(beta.init,lambda.best,folds,datosexplicativos,respuesta,ntau=100){
  
  ngrupos=length(folds)  
  conjtau=seq(0,1,length.out = ntau)  
  mse=rep(0,ntau)
  
  #recorro todos los grupos
  for(k in 1:ngrupos){
    
    datostrain=datosexplicativos[-folds[[k]],]
    datostest=datosexplicativos[folds[[k]],]  
    respuestatrain=respuesta[-folds[[k]]]
    respuestatest=respuesta[folds[[k]]]
    
    ntrain=nrow(datostrain)
    ntest=nrow(datostest)
    
    #recorro todos los folds
    for(j in 1:ntau){
      
      tau=conjtau[j]
      coefthreslassorefit=cv.threslasso(beta.init,datostrain,respuestatrain,tau)$coefthres      
      predichos=predecir(datostest,coefthreslassorefit)
      
      mse[j]=mse[j]+mean((predichos-respuestatest)^2)
    }
    
  }
  tauganador=conjtau[which.min(mse)]
  
  #calculo coeficientes con el tau ganador
  coefthreslassorefitgana=threslasso(beta.init,datostrain,respuestatrain,tauganador)$coefthres
  
  return(list(coef=coefthreslassorefitgana,tau=tauganador,mse=mse))
}








