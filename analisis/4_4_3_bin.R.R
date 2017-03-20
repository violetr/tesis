
packages=c('isingFit','IsingSampler','gRapHD')

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}


source('./metodos/seleccion_modelo_binario.R')
source('./auxiliares.R')

#############################################################

#auxiliares especificos para este dataset

convertirAVariable<-function(i,j,p){return(p*(i-1)+j)}

#crea matriz de grafo estrella con grado d. requiere 0<d<p
crearMatrizEstrella<-function(p,d,w,mixed=FALSE){
  MatrizEstrella=matrix(rep(0,p*p),p)
  vecinosde1=sample(2:p,d)
  for(i in 1:d){
    MatrizEstrella[1,vecinosde1[i]]=w*(sample(c(1,-1),1))^mixed
  }
  return(MatrizEstrella+t(MatrizEstrella))
}

crearMatrizGrilla4n<-function(p,w,mixed=FALSE){
  MatrizGrilla=matrix(rep(0,p^2*p^2),p^2)
  for(i in 1:p){
    for(j in 1:p){
      if(j != i){
        if(i !=1){
          MatrizGrilla[convertirAVariable(i,j,p),convertirAVariable(i-1,j,p)]=(sample(c(1,-1),1)^(mixed))*w
        }
        if(i !=p){
          MatrizGrilla[convertirAVariable(i,j,p),convertirAVariable(i+1,j,p)]=(sample(c(1,-1),1)^(mixed))*w
        }
        if(j !=1){
          MatrizGrilla[convertirAVariable(i,j,p),convertirAVariable(i,j-1,p)]=(sample(c(1,-1),1)^(mixed))*w
        }
        if(j !=p){
          MatrizGrilla[convertirAVariable(i,j,p),convertirAVariable(i,j+1,p)]=(sample(c(1,-1),1)^(mixed))*w
        }
      }
    }
  }
  return(MatrizGrilla)
}


#############################################################

# simu grilla #
 
p=8
w=0.5
thres=rep(0.5,p)
nrep=50
ncasostestigo=5
casostestigo=sort(sample(1:nrep,ncasostestigo))
adycasostestigo=list()
MM=crearMatrizGrilla4n(p,w)
matrizgrilla=MM+t(MM)

adytrue=matrizgrilla !=0 
for(i in 1:dim(adytrue)[1]){
  adytrue[i,i]=FALSE
} 

pdf('./resultados/grilla4n.pdf')
grilla<- network::network(adytrue, directed = FALSE)
ggnet2(grilla,edge.color = "black",color="black",mode="kamadakawai")
dev.off()

n=15000
listamuestrasgrilla=list()
for(i in 1:nrep){
  print(i)
  listamuestrasgrilla[[i]]=IsingSampler(n,matrizgrilla,thres,method = "CFTP")
}

aciertos=0
l=1
for(j in 1:nrep){
  print(j)  
  muestra=listamuestrasgrilla[[j]]
  variablesselec=apply(muestra,2,sum)>5 & nrow(muestra)-apply(muestra,2,sum)>5
  estimoEBIC=IsingFit(muestra[,variablesselec],gamma=0,plot=FALSE)
  ady=estimoEBIC$weiadj !=0 
  for(i in 1:dim(ady)[1]){
    ady[i,i]=FALSE
  }
  
  if(j %in% casostestigo){
    adycasostestigo[[l]]=ady
    l=l+1
  }
  
  newp=nrow(ady)
  estimobien=(newp*newp-sum(ady==adytrue[variablesselec,variablesselec])==0)
  
  aciertos=aciertos+as.integer(estimobien)
}

prueba<- network::network(adycasostestigo[[3]], directed = FALSE)
ggnet2(prueba,edge.color = "black",color="grey",mode="kamadakawai")


# simu estrella #

p=20
d=4
w=0.5
conjn=500
conjn=c(500,1000,5000,10000,15000,20000)
nrep=50
matrizmuestra=crearMatrizEstrella(p,d,w)
thres=rep(0.5,p)

adytrue=matrizmuestra !=0 
for(i in 1:dim(adytrue)[1]){
  adytrue[i,i]=FALSE
} 

estrella<- network::network(matrizmuestra, directed = FALSE)
ggnet2(estrella,edge.color = "black",label=TRUE)

aciertos=matrix(rep(0,2*length(conjn)),nrow=2)

n=50

for(k in 1:length(conjn)){
  
  n=conjn[k]
  print(n)
  
  for(j in 1:nrep){
    
    muestra=IsingSampler(n,matrizmuestra,thres)
    datosfact=apply(muestra,2,as.factor)
    chowliu=minForest(data.frame(datosfact),homog=TRUE,forbEdges=NULL,stat="BIC")
    m=attributes(chowliu)$edges
    attr(m,"n") = ncol(muestra)
    estimacion<- network::network(m, matrix.type="edgelist", directed = FALSE)
    ggnet2(estimacion,label=TRUE)
    matrizestimada<-as.matrix.network.adjacency(estimacion)
    
    adycl=matrizestimada!=0
    for(i in 1:dim(adycl)[1]){
      adycl[i,i]=FALSE
    }
    
    estimobienCL=(p*p-sum(adycl==adytrue)==0)
    
    variablesselec=apply(muestra,2,sum)>5 & nrow(muestra)-apply(muestra,2,sum)>5
    estimoEBIC=IsingFit(muestra[,variablesselec],gamma=0.25,plot=FALSE)
    adyNLR=estimoEBIC$weiadj !=0 
    for(i in 1:dim(ady)[1]){
      adyNLR[i,i]=FALSE
    }
    
    newp=nrow(adyNLR)
    estimobienNLR=(newp*newp-sum(adyNLR==adytrue[variablesselec,variablesselec])==0)
    
    aciertos[1,k]=aciertos[1,k]+as.integer(estimobienCL)
    aciertos[2,k]=aciertos[2,k]+as.integer(estimobienNLR)
    
    variablesselec=apply(muestra,2,sum)>5 & nrow(muestra)-apply(muestra,2,sum)>5
    adyCV=nodewiselogreg(muestra[,variablesselec],and=FALSE,metodo="EBIC",gamma=0.5)
    for(i in 1:dim(adyCV)[1]){
      adyCV[i,i]=FALSE
    }
    
    ver<-network::network(adyCV)
    ggnet2(ver,label=TRUE)
    
    newp=nrow(adyCV)
    estimobienNLRCV=(newp*newp-sum(adyCV==adytrue[variablesselec,variablesselec])==0)
    
    aciertosCV[k]=aciertosCV[k]+as.integer(estimobienNLRCV)
  }
}
aciertos/50

estrella<- network::network(attributes(chowliu)$edges, directed = FALSE)
ggnet2(estrella,edge.color = "lightgrey",label=TRUE)




############################# GRAFICO ###############################

pdf('./resultados/simudiscretas1.pdf')

aciertosNLR<-c(0,0, 0.02, 0.98, 0.96,  1.0, 1.00)
aciertosCL<-c(0,0.02,0.08,0.52,0.52,0.7,0.76)

plot(conjn,aciertosNLR,col="blue",type="b",xlab="n",ylab="Probabilidad de acierto")
points(conjn,aciertosCL,col="red",type="b")

dev.off()

#####################################################################








