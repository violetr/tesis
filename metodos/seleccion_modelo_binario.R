
library(glmnet)

source('./auxiliares.R')

# estimacion modelo binario

# auxiliar

EBIC<-function(n,p,logl,gamma,tamJ){
  -2*logl+tamJ*log(n)+2*gamma*tamJ*log(p-1)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# FUNCION
# nodewiselogreg(datos,and=FALSE,metodo="CV")

# estima el grafo del modelo binario (ising) con el 
# metodo nodewise logistic regression que consiste en 
# realizar regresiones log\'isticas para estimar
# el vecindario de cada nodo

# INPUT:
# datos: Matriz de datos dim=nxp
# and: Regla para establecer el vecindario (AND o OR)
# metodo: metodo para la seleccion del parametro
#         de penalizacion lambda. opciones:
#         c("CV","EBIC")

# OUTPUT:
# matriz de parametros estimados Theta
# matriz con estructura del grafo con booleans

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


#muestras=list()
#muestras[[i]]
#nlambda=50

nodewiselogreg<-function(datos,and=FALSE,metodo="CV",gamma=0.25,nlambda=100){
  n=nrow(datos)
  p=ncol(datos)
  
  coeficientes=matrix(rep(1,p*p),p,p)
  
  for(i in 1:p){
    print(i)
    datosexplicativos=datos[,-i]
    if(metodo=="CV"){
      model=cv.glmnet(datosexplicativos,datos[,i],family="binomial",nfolds=5)
      coeff=coef(model,s="lambda.min")
    }
    if(metodo=="EBIC"){
      model=glmnet(datosexplicativos,datos[,i],family="binomial",nlambda = nlambda)
      nlambda=length(model$lambda)
      EBICm=rep(EBIC(n,p,0,0.25,p-1),nlambda)
      
      for(kiter in 1:nlambda){
        print(kiter)
        variablesseleccionadas=obtenerindices(model$beta[,kiter]!=0)
        if((length(variablesseleccionadas)>1) | (length(variablesseleccionadas)==1 & variablesseleccionadas[1]!=0)){
          
          subdatos=as.data.frame(datosexplicativos[,variablesseleccionadas])
          modellogl=glm(datos[,i] ~ ., family="binomial", data=subdatos)
          logl=logLik(modellogl)
          EBICm[kiter]=EBIC(n,p,logl,gamma,length(variablesseleccionadas))
        }
        
      }
      
      lambdaselec=model$lambda[which.min(EBICm)]
      
      coeff=coef(model,s=lambdaselec)
      
    }
    if(metodo=="stability"){
      lambda=17*sqrt(log(p)/n)
    }
    
    for(j in 1:p){
      #esto lo tengo que separar porque en el vector de coeficientes tengo un lugar menos
      if(j<i){
        coeficientes[i,j]=coeff[j+1]
      }
      if(i<j){
        coeficientes[i,j]=coeff[j]
      }
      if(i==j){
        coeficientes[i,j]=coeff[1]
      }
    }
  }
  AMi=coeficientes !=0
  
  #los de la diagonal quedan en TRUE
  AM=AMi
  for(i in 1:p){
    for(j in 1:p){
      if(i != j){
        if(and==TRUE){
          AM[i,j]=AMi[i,j] && AMi[j,i]
        }else{
          AM[i,j]=AMi[i,j] || AMi[j,i]
        }
      }
    }
  }
  
  cnoceros=log(sum(AM==TRUE))
  
  return(AM)
}
