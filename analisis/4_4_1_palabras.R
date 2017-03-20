#########################bibliotecas#########################

library(R.matlab)
library(glmnet)
library(IsingFit)
library(igraph)

source('C:/Users/violeta/Dropbox/codR/auxiliares.R')

#############################################################

#auxiliares especificos para este dataset

# obtiene la palabla correspondiente 
damePalabra<-function(wordlist,nrosvar,indice){
  return(news$wordlist[[nrosvar[indice]]][[1]][1,1])
}

dameListaPalabras<-function(wordlist,nrosvar){
  lista=rep("",length(nrosvar))
  for( i in 1:length(lista)){
    lista[i]=damePalabra(wordlist,nrosvar,i)
  }
  return(lista)
}

EBIC<-function(n,p,logl,gamma,tamJ){
  -2*logl+tamJ*log(n)+2*gamma*tamJ*log(p-1)
}

############################datos##############################

news<-readMat("C:\\Users\\violeta\\Downloads\\20news_w100.mat",sparseMatrixClass="matrix")

summary(news)

# en documents esta la matrz de datos
# en wordlist la lista de palabras (100)

set.seed(2017)

#me quedo con solo 1000 textos
nrosdoc=sample(dim(news$documents)[2],1000)
parte=news$documents[,nrosdoc]
#transpongo
traspp=t(parte)

# me quedo solo con las palabras que tengan mas de 8 apariciones
# glmnet devuelve warning en caso contrario y creo falla si no hay ninguna
# para kfold - cros validation
nrosvar=obtenerindices(apply(traspp,2,sum)>8)
traspp=traspp[,apply(traspp,2,sum)>8]
dim(traspp)

# vuelvo binarios los bools
datos=traspp*1

nombrescolumnas=dameListaPalabras(datos,nrosvar)
colnames(datos)=nombrescolumnas

#############################################################

modelo=IsingFit(datos,gamma=0.25)

ady=modelo$weiadj !=0 
for(i in 1:dim(ady)[1]){
  ady[i,i]=FALSE
}

# borro los aislados
sinisolated2=ady[apply(ady,2,sum)!=0,apply(ady,2,sum)!=0]
dim(sinisolated2)

net4<- network::network(sinisolated2, directed = FALSE)
ggnet2(net4, size=8, label = TRUE, alpha = 1,label.size =3,color = "lightcyan",edge.color = "lightgrey")


# grafico
net3 <- graph.adjacency(sinisolated2,mode="undirected")
l <- layout.fruchterman.reingold(net3)
l <-norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net3,rescale=FALSE,vertex.color= 'aquamarine2',vertex.frame.color="grey",layout=l*1.5,asp=0.7,vertex.label.color="black",vertex.label.cex=0.8)

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
#         c("CV","wain","EBIC")

# OUTPUT:
# matriz de parametros estimados Theta
# matriz con estructura del grafo con booleans

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


muestras=list()


muestras[[i]]
nlambda=50
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


AM=nodewiselogreg(datos,TRUE,"EBIC",gamma=0.00001,nlambda=50)
dim(AM)
colnames(AM)=nombrescolumnas
rownames(AM)=nombrescolumnas

for(i in 1:dim(AM)[1]){
  AM[i,i]=FALSE
}

sinisolated=AM[apply(AM,2,sum)!=0,apply(AM,2,sum)!=0]

dim(sinisolated)


holi=apply(subgrafo,2,sum)!=0
subgrafo=subgrafo[holi,holi]

net2<- network::network(sinisolated, directed = FALSE)

#net2 <- graph.adjacency(sinisolated,mode="undirected")
#net2 <- graph.adjacency(subgrafo,mode="undirected")
nombrescolumnas

ggnet2(net2, size=8, label = TRUE, alpha = 1,label.size =3,color = "lightcyan",edge.color = "lightgrey")


plot(net2,vertex.color= 'white',vertex.label.color= "black" , vertex.label.dist=0,vertex.label.cex=0.5,edge.curved=0,vertex.size=15,     edge.arrow.size=0.2,
     edge.color="grey",
     edge.width=1,layout= layout.lgl )

l <- layout.fruchterman.reingold(net2)
l <-norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net2,rescale=FALSE,vertex.color= 'aquamarine2',vertex.frame.color="grey",layout=l*1.5,vertex.label.color="black",vertex.label.cex=0.8,asp=0.69)




