
packages=c('IsingFit','IsingSampler','gRapHD')

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

convertirAVariable <- function(i, j, p) return(p * (i - 1) + j)

#crea matriz de grafo estrella con grado d. requiere 0<d<p
crearMatrizEstrella <- function(p, d, w, mixed = FALSE) {
  MatrizEstrella = matrix(rep(0, p * p), p)
  vecinosde1 = sample(2:p, d)
  for (i in 1:d) {
    MatrizEstrella[1, vecinosde1[i]] = w * (sample(c(1, -1), 1))^mixed
  }
  return(list(matriz=MatrizEstrella + t(MatrizEstrella),vecinosde1=vecinosde1))
}

crearMatrizGrilla4n <- function(p, w, mixed = FALSE) {
  MatrizGrilla = matrix(rep(0, p^2 * p^2), p^2)
  for (i in 1:p) {
    for (j in 1:p) {
      if (j != i) {
        if (i != 1) {
          MatrizGrilla[convertirAVariable(i, j, p), convertirAVariable(i - 1, j, p)] = (sample(c(1, -1), 1)^(mixed)) * w
        }
        if (i != p) {
          MatrizGrilla[convertirAVariable(i, j, p), convertirAVariable(i + 1, j, p)] = (sample(c(1, -1), 1)^(mixed)) * w
        }
        if (j != 1) {
          MatrizGrilla[convertirAVariable(i, j, p), convertirAVariable(i, j - 1, p)] = (sample(c(1, -1), 1)^(mixed)) * w
        }
        if (j != p) {
          MatrizGrilla[convertirAVariable(i, j, p), convertirAVariable(i, j + 1, p)] = (sample(c(1, -1), 1)^(mixed)) * w
        }
      }
    }
  }
  return(MatrizGrilla)
}

CalcularPrecisionyRecall <- function(matriztrue,matrizest){
  TP = sum((matrizest==matriztrue) * (matriztrue==TRUE))
  FP = sum((matrizest!=matriztrue) * (matriztrue==FALSE))
  TN = sum((matrizest==matriztrue) * (matriztrue==FALSE))
  FN = sum((matrizest!=matriztrue) * (matriztrue==TRUE))
  
  Precision = TP/(TP+FP)
  Recall = TP/(TP+FN)
  return(list(Precision=Precision,Recall=Recall))  
}


#############################################################

# simu grilla #
 
p = 8
w=0.5
thres = rep(0.5, p)

MM=crearMatrizGrilla4n(p, w)
matrizgrilla = MM + t(MM)

adytrue = matrizgrilla != 0 
for (i in 1:dim(adytrue)[1]) {
  adytrue[i, i] = FALSE
} 


conjn=c(500,1000,5000,10000)

set.seed(2009)

adytrue=matrizgrilla !=0 
for(i in 1:dim(adytrue)[1]){
  adytrue[i,i]=FALSE
} 

nrep=20

n=500
listamuestrasgrilla1=list()
for(i in 1:nrep){
  print(i)
  listamuestrasgrilla1[[i]]=IsingSampler(n,matrizgrilla,thres)
}

n=1000
listamuestrasgrilla2=list()
for(i in 1:nrep){
  print(i)
  listamuestrasgrilla2[[i]]=IsingSampler(n,matrizgrilla,thres)
}

n=5000
listamuestrasgrilla3=list()
for(i in 1:nrep){
  print(i)
  listamuestrasgrilla3[[i]]=IsingSampler(n,matrizgrilla,thres)
}

n=10000
listamuestrasgrilla4=list()
for(i in 1:nrep){
  print(i)
  listamuestrasgrilla4[[i]]=IsingSampler(n,matrizgrilla,thres)
}

grilla<- network::network(matrizgrilla, directed = FALSE)
ggnet2(grilla,edge.color = "black",label=TRUE)

aciertos=matrix(rep(0,4*length(conjn)),nrow=4)
Precision=matrix(rep(0,4*length(conjn)),nrow=4)
Recall=matrix(rep(0,4*length(conjn)),nrow=4)

ncasostestigo = 3
casostestigo = sort(sample(1:nrep, ncasostestigo))
adycasostestigo = list(list(list(),list(),list(),list()),list(list(),list(),list(),list()),list(list(),list(),list(),list()),list(list(),list(),list(),list()))

for(k in 1:length(conjn)){
  
  n=conjn[k]
  print(n)

  nrotestigo=1
  for(j in 1:nrep){
    
    print(paste(n,j,sep=","))
    
    if(n==500){
      muestra=listamuestrasgrilla1[[j]]  
    }
    
    if(n==1000){
      muestra=listamuestrasgrilla2[[j]]  
    }
    
    if(n==5000){
      muestra=listamuestrasgrilla3[[j]]  
    }
    
    if(n==10000){
      muestra=listamuestrasgrilla4[[j]]  
    }
    
    variablesselec=apply(muestra,2,sum)>15 & nrow(muestra)-apply(muestra,2,sum)>15
    estimoEBIC=IsingFit(muestra[,variablesselec],gamma=0.25,plot=TRUE)
    adyNLR=estimoEBIC$weiadj !=0 
    for(i in 1:dim(adyNLR)[1]){
      adyNLR[i,i]=FALSE
    }
    
    newp=nrow(adyNLR)
    estimobienNLR=(newp*newp-sum(adyNLR==adytrue[variablesselec,variablesselec])==0)
    
    adyCV=nodewiselogreg(muestra[,variablesselec],and=TRUE,metodo="CV",gamma=0.25)
    for(i in 1:dim(adyCV)[1]){
      adyCV[i,i]=FALSE
    }
    
    estimobienCV=(newp*newp-sum(adyCV==adytrue[variablesselec,variablesselec])==0)
    
    adystab=nodewiselogreg(muestra[,variablesselec],and=FALSE,metodo="stability",gamma=0.25)
    for(i in 1:dim(adystab)[1]){
      adystab[i,i]=FALSE
    }
    
    estimobienstab=(newp*newp-sum(adystab==adytrue[variablesselec,variablesselec])==0)
    
    aciertos[2,k]=aciertos[2,k]+as.integer(estimobienNLR)
    aciertos[3,k]=aciertos[3,k]+as.integer(estimobienCV)
    aciertos[4,k]=aciertos[4,k]+as.integer(estimobienstab)
    
    Precision[2,k]=Precision[2,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyNLR)$Precision
    Precision[3,k]=Precision[3,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyCV)$Precision
    Precision[4,k]=Precision[4,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adystab)$Precision
    
    Recall[2,k]=Recall[2,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyNLR)$Recall
    Recall[3,k]=Recall[3,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyCV)$Recall
    Recall[4,k]=Recall[4,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adystab)$Recall
    
    if (j %in% casostestigo) {
      
      adycasostestigo[[2]][[k]][[nrotestigo]]=adyNLR
      adycasostestigo[[3]][[k]][[nrotestigo]]=adyCV
      adycasostestigo[[4]][[k]][[nrotestigo]]=adystab
      nrotestigo=nrotestigo+1
    }
  }
}


# hidden = matrix(rep(0, 4 * 64 * 64), 64 * 2)
# 
# dim(hidden)
# 
# hidden[1:64, 1:64] = adytrue
# 
# hidden[1:64, 65:128] = diag(64)
# 
# hidden[65:128,1:64] = diag(64)
# 
# hidden[65:128,65:128] = matrix(rep(0, 64*64), 64)
# 
# pdf('./resultados/grilla4n.pdf')
# grilla <- network::network(adytrue, directed = FALSE)
# ggnet2(grilla, size = 2.5, mode = "kamadakawai", color = "black", edge.color = "lightgrey")
# dev.off()
# 
# pdf('./resultados/grilla4nhidden.pdf')
# hidden <- network::network(hidden, directed = FALSE)
# coloresh <- c("black", "darkslategray4")
# names(coloresh) = c("black", "darkslategray4")
# hidden %v% "color" <- c(rep("black", 64), rep("darkslategray4", 64))
# ggnet2(hidden,size = 3,mode = "kamadakawai", color = coloresh[hidden %v% "color"], edge.color = "lightgrey")
# dev.off()
# 
# hidden[ , 65]

# n=500
# listamuestrasgrilla=list()
# for(i in 1:nrep){
#   print(i)
#   listamuestrasgrilla[[i]]=IsingSampler(n,matrizgrilla,thres,method = "CFTP")
# }
# 
# aciertos=0
# l=1
# for(j in 1:nrep){
#   print(j)  
#   muestra=listamuestrasgrilla[[j]]
#   variablesselec=apply(muestra,2,sum)>5 & nrow(muestra)-apply(muestra,2,sum)>5
#   estimoEBIC=IsingFit(muestra[,variablesselec],gamma=0,plot=FALSE)
#   ady=estimoEBIC$weiadj !=0 
#   for(i in 1:dim(ady)[1]){
#     ady[i,i]=FALSE
#   }
#   
#   if(j %in% casostestigo){
#     adycasostestigo[[l]]=ady
#     l=l+1
#   }
#   
#   newp=nrow(ady)
#   estimobien=(newp*newp-sum(ady==adytrue[variablesselec,variablesselec])==0)
#   
#   aciertos=aciertos+as.integer(estimobien)
# }

prueba<- network::network(adycasostestigo[[3]], directed = FALSE)
ggnet2(prueba,edge.color = "black",color="grey",mode="kamadakawai")


# simulacion estrella #

p=64
d=7
w=0.5

conjn=c(500,1000,5000,10000)

set.seed(2009)

CREARESTRELLA <- crearMatrizEstrella(p,d,w)
matrizestrella=CREARESTRELLA$matriz
thres=rep(0.2,p)

vecinosde1 <- CREARESTRELLA$vecinosde1

adytrue=matrizestrella !=0 
for(i in 1:dim(adytrue)[1]){
  adytrue[i,i]=FALSE
} 

nrep=20

n=500
listamuestrasestrella1=list()
for(i in 1:nrep){
  print(i)
  listamuestrasestrella1[[i]]=IsingSampler(n,matrizestrella,thres,method = "CFTP")
}

n=1000
listamuestrasestrella2=list()
for(i in 1:nrep){
  print(i)
  listamuestrasestrella2[[i]]=IsingSampler(n,matrizestrella,thres,method = "CFTP")
}

n=5000
listamuestrasestrella3=list()
for(i in 1:nrep){
  print(i)
  listamuestrasestrella3[[i]]=IsingSampler(n,matrizestrella,thres,method = "CFTP")
}

n=10000
listamuestrasestrella4=list()
for(i in 1:nrep){
  print(i)
  listamuestrasestrella4[[i]]=IsingSampler(n,matrizestrella,thres,method = "CFTP")
}


estrella<- network::network(matrizestrella, directed = FALSE)
ggnet2(estrella,edge.color = "black",label=TRUE)

aciertos=matrix(rep(0,4*length(conjn)),nrow=4)
Precision=matrix(rep(0,4*length(conjn)),nrow=4)
Recall=matrix(rep(0,4*length(conjn)),nrow=4)

ncasostestigo = 3
casostestigo = sort(sample(1:nrep, ncasostestigo))
adycasostestigo = list(list(list(),list(),list(),list()),list(list(),list(),list(),list()),list(list(),list(),list(),list()),list(list(),list(),list(),list()))

for(k in 2:length(conjn)){
  
  n=conjn[k]
  print(n)
  nrotestigo=1
  n=500
  k=1
  for(j in 1:nrep){
    
    print(paste(n,j,sep=","))
    
    if(n==500){
      muestra=listamuestrasestrella1[[j]]  
    }

    if(n==1000){
      muestra=listamuestrasestrella2[[j]]  
    }

    if(n==5000){
      muestra=listamuestrasestrella3[[j]]  
    }

    if(n==10000){
      muestra=listamuestrasestrella4[[j]]  
    }
    
    datosfact=apply(muestra,2,as.factor)
    chowliu=minForest(data.frame(datosfact),homog=TRUE,forbEdges=NULL,stat="BIC")
    m=attributes(chowliu)$edges
    attr(m,"n") = ncol(muestra)
    estimacion<- network::network(m, matrix.type="edgelist", directed = FALSE)
    ggnet2(estimacion,label=TRUE)
    matrizestimada<-as.matrix.network.adjacency(estimacion)
    # 
    adycl=matrizestimada!=0
    for(i in 1:dim(adycl)[1]){
      adycl[i,i]=FALSE
    }
    
    estimobienCL=(p*p-sum(adycl==adytrue)==0)
     
    variablesselec=apply(muestra,2,sum)>5 & nrow(muestra)-apply(muestra,2,sum)>5
    estimoEBIC=IsingFit(muestra[,variablesselec],gamma=0.25,plot=TRUE)
    adyNLR=estimoEBIC$weiadj !=0 
    for(i in 1:dim(adyNLR)[1]){
      adyNLR[i,i]=FALSE
    }

    newp=nrow(adyNLR)
    estimobienNLR=(newp*newp-sum(adyNLR==adytrue[variablesselec,variablesselec])==0)
 
    adyCV=nodewiselogreg(muestra[,variablesselec],and=TRUE,metodo="CV",gamma=0.25)
    for(i in 1:dim(adyCV)[1]){
       adyCV[i,i]=FALSE
    }
    
    estimobienCV=(newp*newp-sum(adyCV==adytrue[variablesselec,variablesselec])==0)
    
    adystab=nodewiselogreg(muestra[,variablesselec],and=FALSE,metodo="stability",gamma=0.25)
    for(i in 1:dim(adystab)[1]){
      adystab[i,i]=FALSE
    }

    estimobienstab=(newp*newp-sum(adystab==adytrue[variablesselec,variablesselec])==0)
  
    aciertos[1,k]=aciertos[1,k]+as.integer(estimobienCL)
    aciertos[2,k]=aciertos[2,k]+as.integer(estimobienNLR)
    aciertos[3,k]=aciertos[3,k]+as.integer(estimobienCV)
    aciertos[4,k]=aciertos[4,k]+as.integer(estimobienstab)
    
    Precision[1,k]=Precision[1,k]+CalcularPrecisionyRecall(adytrue,adycl)$Precision
    Precision[2,k]=Precision[2,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyNLR)$Precision
    
    Precision[3,k]=Precision[3,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyCV)$Precision
    Precision[4,k]=Precision[4,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adystab)$Precision

    Recall[1,k]=Recall[1,k]+CalcularPrecisionyRecall(adytrue,adycl)$Recall
    Recall[2,k]=Recall[2,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyNLR)$Recall
    Recall[3,k]=Recall[3,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adyCV)$Recall
    Recall[4,k]=Recall[4,k]+CalcularPrecisionyRecall(adytrue[variablesselec,variablesselec],adystab)$Recall
    
    if (j %in% casostestigo) {
      
      adycasostestigo[[1]][[k]][[nrotestigo]]=adycl
      adycasostestigo[[2]][[k]][[nrotestigo]]=adyNLR
      adycasostestigo[[3]][[k]][[nrotestigo]]=adyCV
      adycasostestigo[[4]][[k]][[nrotestigo]]=adystab
      nrotestigo=nrotestigo+1
    }
  }
}


Recall/10
Precision/10

estrella<- network::network(attributes(chowliu)$edges, directed = FALSE)
ggnet2(estrella,edge.color = "lightgrey",label=TRUE)

estimacion<- network::network(adycasostestigo[[4]][[2]][[1]], directed = FALSE)
ggnet2(estimacion,label=TRUE)

############################# GRAFICO ###############################

pdf('./resultados/simudiscretas1.pdf')

aciertosNLR<-c(0,0, 0.02, 0.98, 0.96,  1.0, 1.00)
aciertosCL<-c(0,0.02,0.08,0.52,0.52,0.7,0.76)

plot(conjn,aciertosNLR,col="blue",type="b",xlab="n",ylab="Probabilidad de acierto")
points(conjn,aciertosCL,col="red",type="b")

dev.off()

#####################################################################








