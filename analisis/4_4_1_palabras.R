#########################bibliotecas#########################

packages=c('R.matlab','IsingFit')

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

############################datos##############################

news<-readMat("./datos/20news_w100.mat",sparseMatrixClass="matrix")

summary(news)

# en documents esta la matrz de datos
# en wordlist la lista de palabras (100)

set.seed(2017)

#me quedo con solo 1000 textos
nrosdoc=sample(dim(news$documents)[2],1000)
nrosdoc2=sample(dim(news$documents)[2],6000)
parte=news$documents[,nrosdoc2]
#transpongo
traspp=t(parte)

# me quedo solo con las palabras que tengan mas de 8 apariciones
# glmnet devuelve warning en caso contrario y creo falla si no hay ninguna
# para kfold - cros validation
nrosvar=obtenerindices(apply(traspp,2,sum)>70)
traspp=traspp[,apply(traspp,2,sum)>70]
dim(traspp)

# vuelvo binarios los bools
datos=traspp*1

nombrescolumnas=dameListaPalabras(datos,nrosvar)
colnames(datos)=nombrescolumnas

#############################################################

modelo=IsingFit(datos,gamma=0.0001)

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


AM2=nodewiselogreg(datos,TRUE,"stability",gamma=0.5,nlambda=50)
AM2=nodewiselogreg(datos,TRUE,"EBIC",gamma=0.5,nlambda=50)
AM=nodewiselogreg(datos,TRUE,"CV",gamma=0.00001,nlambda=50)
dim(AM2)
dim(AM)
colnames(AM)=nombrescolumnas
rownames(AM)=nombrescolumnas

for(i in 1:dim(AM)[1]){
  AM[i,i]=FALSE
}

sinisolated=AM[apply(AM,2,sum)!=0,apply(AM,2,sum)!=0]

dim(sinisolated)

colnames(AM2)=nombrescolumnas
rownames(AM2)=nombrescolumnas

for(i in 1:dim(AM2)[1]){
  AM2[i,i]=FALSE
}

sinisolated2=AM2[apply(AM2,2,sum)!=0,apply(AM2,2,sum)!=0]

dim(sinisolated2)


holi=apply(subgrafo,2,sum)!=0
subgrafo=subgrafo[holi,holi]

net2<- network::network(sinisolated, directed = FALSE)

#net2 <- graph.adjacency(sinisolated,mode="undirected")
#net2 <- graph.adjacency(subgrafo,mode="undirected")
nombrescolumnas

ggnet2(net2, size=8, label = TRUE, alpha = 1,label.size =3,color = "lightcyan",edge.color = "lightgrey",mode="kamadakawai")

net22<- network::network(sinisolated2, directed = FALSE)

#net2 <- graph.adjacency(sinisolated,mode="undirected")
#net2 <- graph.adjacency(subgrafo,mode="undirected")
nombrescolumnas

ggnet2(net22, size=8, label = TRUE, alpha = 1,label.size =3,color = "lightcyan",edge.color = "lightgrey",mode="kamadakawai")

plot(net2,vertex.color= 'white',vertex.label.color= "black" , vertex.label.dist=0,vertex.label.cex=0.5,edge.curved=0,vertex.size=15,     edge.arrow.size=0.2,
     edge.color="grey",
     edge.width=1,layout= layout.lgl )

l <- layout.fruchterman.reingold(net2)
l <-norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net2,rescale=FALSE,vertex.color= 'aquamarine2',vertex.frame.color="grey",layout=l*1.5,vertex.label.color="black",vertex.label.cex=0.8,asp=0.69)




