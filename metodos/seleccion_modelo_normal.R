# METODOS seleccion modelo normal

# estimacion del grafo normal:
# Adaptive GLasso: adaptive.glasso
# (Adaptive)GLasso CV: cv.glasso
# Nodewise regression: 
# Stability Selection:
# Gelato (nodewise+threshold):


# metodo adaptive glasso
adaptive.glasso<-function(strain,lambda,penalizardiagonal = FALSE){
  estimaglasso<-glasso(strain,rho=lambda,penalize.diagonal = penalizardiagonal)$wi
  ceros<-indicescero(estimaglasso!=0)
  if(is.null(dim(ceros))){
    matrizpenalizacion<-lambda*ifelse(is.infinite(1/abs(estimaglasso)),0,1/abs(estimaglasso))
    diag(matrizpenalizacion)<-rep(0,nrow(matrizpenalizacion))
    estimaadaptive<-glasso(strain,rho=matrizpenalizacion,penalize.diagonal = penalizardiagonal)$wi     
  }else{
    matrizpenalizacion<-lambda*ifelse(is.infinite(1/abs(estimaglasso)),0,1/abs(estimaglasso))
    diag(matrizpenalizacion)<-rep(0,nrow(matrizpenalizacion))
    estimaadaptive<-glasso(strain,rho=matrizpenalizacion,zero=ceros,penalize.diagonal = penalizardiagonal)$wi
  }
  return(estimaadaptive)
}


#glasso cross-validation
#selecciono lambda glasso

cv.glasso<-function(datos,grupos,conjlambda=seq(0.000001,1,0.005),dibujar=FALSE,adaptive=FALSE){
  
  n=nrow(datos)
  k=length(grupos)
  
  neggauss<-rep(0,length(conjlambda))
  lognoceros<-rep(0,length(conjlambda))
  
  for (i in 1:k){
    
    nrostest<-grupos[[i]]
    
    valid<-datos[nrostest,]
    train<-datos[-nrostest,]
    
    svalid<-estimadors(valid)
    strain<-estimadors(train)
    
    for (j in 1:length(neggauss)){
      
      lambda=conjlambda[j]
      if(adaptive){
        estimaglasso<-adaptive.glasso(strain,lambda,penalizardiagonal = FALSE)  
      }else{
        estimaglasso<-glasso(strain,rho=lambda,penalize.diagonal = FALSE)$wi  
      }
      neggauss[j]<-neggauss[j]+negativegaussloss(estimaglasso,valid,lambda)
      lognoceros[j]<-log(sum(estimaglasso!=0))
    }
    
  }
  
  lambda.gana=conjlambda[which.min(neggauss)]
  
  s<-estimadors(datos)
  
  patronceros=grafo(s,lambda.gana,dibujar)
  
  return(list(lambda.gana,patronceros,neggauss,lognoceros))
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# FUNCION
# nodewisereg(datos,and=FALSE,metodo="CV"):
# estima el grafo del modelo gaussiano con el 
# metodo nodewise (thresholded) regression que consiste en 
# realizar regresiones para estimar
# el vecindario de cada nodo (Meinshausen & Buhlmann 2006)

# INPUT:
# datos: Matriz de datos dim=nxp
# and: Regla para establecer el vecindario (AND o OR)
# metodo: metodo para la seleccion del parametro
#         de penalizacion lambda. opciones:
#         c("CV","wain","EBIC")

# OUTPUT:
# matriz de parametros estimados Theta
# matriz con estructura del grafo con booleans
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
nodewisereg<-function(datos,tau,and=TRUE,lambda,thresholding=TRUE){
  n=nrow(datos)
  p=ncol(datos)
  
  coeficientes=matrix(rep(1,p*p),p,p)
  
  for(i in 1:p){
    
    datosexplicativos=datos[,-i]
    if(thresholding){
      coef=threslasso(datosexplicativos,respuesta=datos[,i],tau,lambda)$coefthres
    }else{
      modelolasso<-glmnet(x=datosexplicativos,y=respuesta,alpha=1,lambda = lambda)
      coef=as.matrix(modelolasso$beta)[,1]
    }
    
    for(j in 1:p){
      #esto lo tengo que separar porque en el vector de coeficientes tengo un lugar menos
      if(j<i){
        coeficientes[i,j]=coef[j]
      }
      if(i<j){
        coeficientes[i,j]=coef[j-1]
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


###############################################################################
######################stability selection graphical models#####################
###############################################################################

#INPUT:
#datos:datos sin estandarizar
#lambda a utilizar, eventualmente cv
#nsubmuestra: tamano de la submuestra de subsampling
#esperanza falsos positivos
#pithr: valor de pi de corte
#OUTPUT:
#grafo estimado

stability.selection.grafico<-function(datos,lambda,B,nsubmuestra=floor(nrow(datos)/2),expfalsepos,pithr,dibujar=TRUE){
  
  n=nrow(datos)
  p=ncol(datos)
  
  #matriz donde guardo los valores de probabilidad estimada
  pi=matrix(rep(0,p*p),ncol=p)
  
  #numero de aristas que voy a tomar como máximo:
  #si ninguna de las elegidas 
  q=floor(sqrt(expfalsepos*p*(2*pithr-1))) 
  
  #voy a tomar B submuestras de tamanio n/2, sin reposicion
  for(k in 1:B){
    
    #sorteo la submuestra
    submuestra=sample(1:n,nsubmuestra,replace=FALSE)
    
    #tomo el valor absoluto de la matriz de precision estimada con la submuestra con glasso
    absprecision<-abs(glasso(var(datos[submuestra,]),lambda,penalize.diagonal = FALSE)$wi)
    
    #pongo cero en la diagonal para no contar a estos lugares como seleccionados
    absprecision[ row(absprecision) == col(absprecision) ]=0
    
    #me quedo con con las q que tengan los valores absolutos de coeficientes mas altos,
    #seleccion es una matrix de p*p con un 1 en las que se estimaron distinto de 0
    #OBS: todas las aristas dentro de las q seleccionadas que tengan valor 0, son desechadas
    #OBS2:tomo las primeras 2*q porque al ser simetrica la matriz si selecciono ij tambien ji
    seleccion<-matrix(as.integer(absprecision %in% head(setdiff(sort(absprecision, TRUE),0), 2*q)),nr = nrow(absprecision))
    
    #para cada lugar de la matriz  sumo uno si ij fue seleccionada
    pi=pi+seleccion
    
  }
  #divido por la cantidad para tener la frecuencia relativa
  pi=pi/B
  
  #indices en una matriz de 2xcantidad de los lugares de la matriz pi mas grandes que el limite pi establecido
  indicesmaxpi <- which(matrix(pi>pithr,nr = nrow(pi)), arr.ind = TRUE)
  
  #vale uno si pi esmayor ue pithr
  nocerosstability<-matrix(rep(0,p*p),ncol=p)
  nocerosstability[indicesmaxpi]=1
  
  #grafo resultante:
  pal3 <- brewer.pal(10, "Set3")
  if(dibujar){
    grafo <- graph.adjacency(nocerosstability,mode="undirected") #uso AM como matriz de adyacencia
    plot(grafo,layout=layout.circle(grafo),vertex.label=NA,vertex.size=6, edge.color="black",vertex.color=rev(pal3))
  }
  
  return(list(pi,nocerosstability))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# FUNCION
# Gelato(datos,grupos,conjlambda=seq(0.000001,1,0.005),tau,and=TRUE):
# estima la inversa de la matriz de covarianza 
# del modelo gaussiano con el 
# metodo gelato (Zhou et al 2011) 
# que consiste en estimar el grafo con
# nodewise thresholded regression y luego estimar por MV
# sin restricciones
# eligiendo el valor de lambda con CV

# INPUT:
# datos: Matriz de datos estandarizados dim=nxp
# grupos: lista con conjuntos(vectores) de indices 
#         para cada fold de CV
# and: Regla para establecer el vecindario (AND o OR)
# conjlambda: vector de valores de lambda a seleccionar
# tau: valor de thresholding
# dibujar: dibujar el grafo estimado

# OUTPUT:
# gana: lambda seleccionado con CV 
# loss: (neg log loglikelihood) de kfolds CV
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
cv.Gelato<-function(datos,grupos,conjlambda=seq(0.000001,1,0.005),tau,and=TRUE){
  
  n=nrow(datos)
  
  k=length(grupos)
  neggauss<-rep(0,length(conjlambda))
  
  for (i in 1:length(grupos)){
    
    nrostest<-grupos[[i]]
    valid<-datos[nrostest,]
    train<-datos[-nrostest,]
    
    for (j in 1:length(neggauss)){
      
      lambda=conjlambda[j]
      
      #aplico nodewise regression a los datos, asi obtengo un grafo estimado
      Matrizbool=nodewisereg(train,tau,and,lambda)
      
      #obtengo los indices que tienen que ser 0 en la matriz de precision
      Indicescero=indicescero(Matrizbool)
      
      #estimador MV de la matriz de covarianza de los datos
      s=var(train)
      
      #estimo lasso sin penalidad y con restriccion de ciertos valores cero,guardo matriz precision estimada
      matrizestimada=glasso(s,rho=0,zero=Indicescero,penalize.diagonal = FALSE)$wi
      
      neggauss[j]<-neggauss[j]+negativegaussloss(matrizestimada,valid,lambda)
    }
    
  }
  
  return(list(loss=neggauss,gana=conjlambda[which.min(neggauss)]))
}

