######################################################################################################################
######################################################################################################################
# AN 3.7.3

#analisis thaliana con glasso y nodewise regression con thresholded lasso distintos tau
#replicar grafico Fig 13.3 buhlmann, pagina 443

#bibliotecas usadas: library(glmnet),library(glasso),library(igraph),library(RColorBrewer)

if (!require(MASS)) {
  install.packages(MASS)
  library(MASS)
}

source('./metodos/seleccion_modelo_normal.R')
source('./auxiliares.R')

######################################################################################################################
############################################# cargo y centro datos ###################################################

#datos biologicos isoprenoid arabidopsis thaliana, n=118 , p=39
#como estimo la matriz de precision 39*39=1521, n=118

#cargo datos
datosthalia=read.table('./datos/datosthaliana.txt',header = TRUE)
matrizthalia=t(datosthalia[,7:length(datosthalia)])

# centro los datos
matrizthalia=scale(matrizthalia)


#####################################################################################################################
#####################################################################################################################

set.seed(2005)

################################################ ESTIMACION CON GLASSO ##############################################

# k fold cross validation

n=nrow(matrizthalia)

grupos=kfoldgrupos(n,10)

# grafico loss glasso

#lambdasprueba=nodewise.getlambdasequence(matriz)
#ak=c(0.0001,0.01,0.05,0.1,0.3,0.5,1,2,4)
#lambdasprueba2=ak*sqrt(log(p)/n)

length(seq(0.000001,1,0.005))
graficoglasso<-function(datos,grupos,  conjlambda=seq(0.000001,1,0.005),
                        lambdaselecc=c(1,2,3,4,5,6,7,10,15,20,53,80,140,180,190),adaptive=FALSE){
  
  ans=cv.glasso(datos,grupos,conjlambda[lambdaselecc],FALSE,adaptive)
  
  s<-estimadors(datos)
  
  nfig=length(lambdaselecc)
  
  lognnozeros=rep(0,nfig)
  neggaussloss=rep(0,nfig)
  for(i in 1:nfig){
    lognnozeros[i]=(ans[[4]])[i]
    neggaussloss[i]=(ans[[3]])[i]
  }
  print(neggaussloss)
  return(list(lognnozeros=lognnozeros,neggaussloss=neggaussloss))
}

g1<-graficoglasso(matrizthalia,grupos,adaptive=TRUE)
plot(g1$lognnozeros,g1$neggaussloss,type="b")
g2<-graficoglasso(matrizthalia,grupos,adaptive=FALSE)
plot(g2$lognnozeros,g2$neggaussloss,type="b",col="red")

################################################# NODEWISE REGRESION ################################################

#grafico loss con nodewise regression

graficonodewisereg<-function(datos,grupos,conjlambda=seq(0.000001,1,0.005),
                             lambdaselecc2=c(1,2,3,10,20,30,40,50,60,80,90,100,120,140,180),
                             conjtau=c(0.012,0.025,0.045,0.062,0.125),
                             and=TRUE){
  
  conjlambda2=conjlambda[lambdaselecc2]

  nfig=length(conjlambda2)
  ntau=length(conjtau)
  
  lognnozerostau=matrix(rep(0,ntau*nfig),nrow=ntau,ncol=nfig)
  neggausslosstau=matrix(rep(0,ntau*nfig),nrow=ntau,ncol=nfig)
  
  for(i in 1:length(conjtau)){
    
    tau=conjtau[i]
    
    node=cv.Gelato(datos,grupos,conjlambda2,tau,and)
    
    for(j in 1:nfig){
      lambda=conjlambda2[j]
      Matrizbool=nodewisereg(datos,tau,and,lambda)
      Indicescero=indicescero(Matrizbool)
      s=var(datos)
      matrizestimada=glasso(s,rho=0,zero=Indicescero,penalize.diagonal = FALSE)$wi
      lognnozerostau[i,j]=log(sum(matrizestimada!=0))
      neggausslosstau[i,j]=(node$loss)[j]
    }
  }
  return(list(lognnozeros=lognnozerostau,neggaussloss=neggausslosstau))
}

g2<-graficonodewisereg(matrizthalia,grupos)
conjtau=c(0.012,0.025,0.045,0.062,0.125)
ntau=length(conjtau)
pal3=brewer.pal(6,"Set1")
for(i in 1:ntau){
  points(g2$lognnozeros[i,],g2$neggaussloss[i,],type="b",col=pal3[i])
}

#####################################################################################################################
####################################################### FIGURA 13.3 #################################################


figura13.3<-function(datos,grupos,conjtau=c(0.012,0.025,0.045,0.062,0.125),and=TRUE){
  
  pal3=brewer.pal(6,"Set1")
    
  n=nrow(datos)
  
  #glasso
  
  g1<-graficoglasso(datos,grupos,adaptive=FALSE)
  g2<-graficoglasso(datos,grupos,adaptive =TRUE)
  
  
  if(and==TRUE){
    titulo="Regla 'AND'"
  }else{
    titulo="Regla 'OR'"
  }
  
  plot(g1$lognnozeros,g1$neggaussloss,type="b",xlab="log(nro!=0matrizprecision)",ylab="negative gaussian loss",main=titulo)
  lines(g2$lognnozeros,g2$neggaussloss,type="b",col="black",pch=2)
  
  #nodewise
  
  g2<-graficonodewisereg(datos,grupos,conjtau=conjtau,and=and)
  
  ntau=length(conjtau)

  for(i in 1:ntau){
    points(g2$lognnozeros[i,],g2$neggaussloss[i,],type="b",col=pal3[i])
  }
  
  neggausslosss=0
  for(k in 1:length(grupos)){
    valid=datos[grupos[[k]],]
    train=datos[-grupos[[k]],]
    neggausslosss=neggausslosss+negativegaussloss(estimadors(train),valid,0)   
  }
  
  print(neggausslosss)
    
  #referencias
  
  legend("bottomleft", 
         c("glasso","adaptive.glasso","tau=0.012","tau=0.025","tau=0.045","tau=0.062","tau=0.125"),
         lwd=1, lty=c(1,1), 
         pch=c(1,2,1,1,1,1,1), 
         cex=0.8, pt.cex = 1,
         col=c('black','black',pal3[1:5]))
}

grupos=kfoldgrupos(n,10)

figura13.3(matrizthalia,grupos,and=FALSE)

# SIMU fig 13.3

p=40
n=180

m=matrix(rnorm(p*p),ncol=p)
Theta=m%*%t(m)

r=sample(1:(p*p),700,replace = FALSE)

for(i in 1:length(r)){
  
  fila=((r[i]-1) %/% p) +1
  columna=((r[i]-1) %% p) +1 
  
  if(fila!=columna){
    Theta[fila,columna]=0
    Theta[columna,fila]=0
  }
}

log(sum(Theta!=0))

Sigma=solve(Theta)
Sigma=round(Sigma,10) # por el error de computo inversa

mu=rep(0,p)
muestra=mvrnorm(n, mu, Sigma)
muestra=scale(muestra)

grupos2=kfoldgrupos(n,10)
figura13.3(muestra,grupos2,and=FALSE)



