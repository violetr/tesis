################## AUXILIARES ###########################

# kfoldgrupos: genera grupos para kfold-CV

# predecir: predict lineal, con vector de coeficientes

# traza: calcula la traza de una matriz

# center_scale: centra los datos y/o los reescala 

# negativegaussloss

# obtenerindices(true): devuelve un vector de los indices
# de un vector que valen true

# indicescero:
# obtiene una matriz con los indices de vectorbool que valen false, 
# los devuelve en una matriz de dos columnas

# estimadors:
# matriz de covarianza empirica

# grafo:
# grafica el grafo correspondiente como circulo y devuelve 
# el objeto grafo y el patron de ceros en la matriz de precision 

# EBIC:
# calcula la penalizacion EBIC

################## BIBLIOTECAS ##########################

packages=c('glasso','glmnet','igraph','RColorBrewer','GGally', 'network')

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

################## FUNCIONES ############################

# dado un n tamano de datos entrenamiento y un
# valor de k, sortea los grupos para hacer 
# k-fold cross validation
kfoldgrupos<-function(n,k=10){
  # n= tamanio muestra train
  # k= cantidad de folds 
  
  sizes<-rep(floor(n/k),k)
  b=n-floor(n/k)*k
  i=1
  while(b>0){
    sizes[i]=sizes[i]+1
    i=i+1
    b=b-1
  }
  
  grupos<-list()
  numeros<-1:n
  
  for(i in 1:k){
    grupos[[i]]<-sample(numeros,sizes[i],replace=FALSE)
    numeros<-setdiff(numeros,grupos[[i]])
  }
  
  return(grupos)
  
}

# obtiene un vector con los indices de vectorbool que valen true,
# si no hay ningun TRUE devuelve el valor "0".
# INPUT: vectorbool = vector de valores bool
# OUTPUT: indices= vector de indices que valen true en vectorbool 
obtenerindices<-function(vectorbool){
  indices=rep(0,1)
  if(sum(vectorbool==TRUE)>0){
    indices=rep(0,sum(vectorbool==TRUE))
    j=1
    for(i in 1:length(vectorbool)){
      if(vectorbool[i]==TRUE){
        indices[j]=i
        j=j+1
      }
    }  
  }
  return(indices)
}
#test:
#obtenerindices(c(FALSE,FALSE,FALSE))

#obtiene una matriz con los indices de matrizbool que valen false, los devuelve en una matriz de dos columnas
#VER que pasa si ninguno es cero!
indicescero<-function(matrizbool){
  p=ncol(matrizbool)
  numerodeceros=sum(matrizbool==FALSE)
  if(numerodeceros>0){
    indices=matrix(rep(0,2*numerodeceros),ncol=2,nrow=numerodeceros)
    k=1
    for(i in 1:p){
      for(j in 1:p){
        if(matrizbool[i,j]==FALSE){
          indices[k,1]=i
          indices[k,2]=j
          k=k+1
        }
      }
    }
  }else{
    indices=0
  }
  return(indices)
}

# mismo resultado que predict (lineal)
# dado una matriz de datos test y un vector de 
# coeficientes beta estimados predice el valor
# y para cada fila de la matriz de datos
predecir<-function(datostest,coefbeta){
  #datostest: matriz de dimension w x p con datos nuevos
  #coefbeta: vector de coeficientes, longitud p 
  return(datostest  %*% coefbeta)
}


#funcion centra pero no reescala
center_scale <- function(x) {
  scale(x, scale = FALSE)
}

#calcula la traza de una matriz
traza<-function(A){
  return(sum(diag(A)))
}

#funcion de perdida a utilizar en la validacion cruzada
#theta es la estimacion de la matriz de precision con el conjunto train
#valid es el conjunto de datos para validacion
negativegaussloss<-function(theta,valid,lambda){
  #return(-log(det(theta))+traza(estimadors(valid)%*%theta))
  return(-log(det(theta))+traza(var(valid)%*%theta))
}


#estimador maxima verosimilitud matriz de covarianza
estimadors<-function(datos){
  s=0
  n=nrow(datos)
  for(i in 1:n){
    s=s+(datos[i,]-apply(datos,2,mean))%*%t(datos[i,]-apply(datos,2,mean)) 
  }
  s=s/n
  return(s)
}
#test: estimadors
#sigma=diag(4)+matrix(c(0,3,0,1,0,0,0,0,0,0,0,0,0,0,0,0),4,4)+t(matrix(c(0,3,0,1,0,0,0,0,0,0,0,0,0,0,0,0),4,4))
#spd=sigma%*%t(sigma)
#datos=rmvnorm(100,rep(0,4),spd)
#estimadors(datos)


#grafica el grafo correspondiente como circulo y devuelve el objeto grafo y el patron de ceros en la matriz de precision 
grafo<-function(covarianza,lambda,dibujar=FALSE){
  
  res.lasso<-glasso(covarianza,rho=lambda,penalize.diagonal = FALSE) 
  
  AM <- res.lasso$wi != 0 #true en donde la matriz de covarianza distinta de 0
  
  diag(AM) <- F # no pongo aristas de que involucren solo un vertice para graficar
  
  grafo <- graph.adjacency(AM,mode="undirected") #uso AM como matriz de adyacencia
  pal3 <- brewer.pal(10, "Set3")
  if(dibujar==TRUE){
    plot(grafo,layout=layout.circle(grafo),vertex.size=8,vertex.label=NA,edge.color="black",vertex.color=rev(pal3))
    
    #plot(grafo,layout=layout.circle(grafo),vertex.label.cex=0.2,vertex.size=15,vertex.label.color= "black" , edge.color="black",vertex.color=rev(pal3))
  }
  return(list(grafo=grafo,matriz=AM))
  
}



