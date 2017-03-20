#######################################################################
##############################bibliotecas #############################

#en esta esta implementado stability selection para regresion, y estan los datos riboflavin
library(hdi)
source('C:/Users/violeta/Dropbox/tesis violeta/codR/glasso.R')

#######################################################################
############################### datos #################################

#carga los datos ubicados en el paquete hdi:
data("riboflavin")
#estan en 2 columnas, respuesta y explicativas

#despliego las explicativas
dataribo<-matrix(rep(0,71*4088),nrow=71,ncol=4088)
for(i in 1:71){
  dataribo[i,]=riboflavin[i,2]
}

#tomo una muestra de 160 variables de las 4088 originales
set.seed(2010)
sorteo=sample(1:4088,160,replace=FALSE)

submuestraribo<-dataribo[,sorteo]

#grupos=kfoldgrupos(71,5)
#resultglassocv=glasso.cv(submuestraribo,grupos,conjlambda=seq(0.051,0.95,0.05),dibujar=FALSE)
 
#me genero los datos con permutaciones independientes, queda todo con covarianza 0
todopermutado<-matrix(rep(0,71*160),nrow=71,ncol=160)
for(j in 1:ncol(submuestraribo)){
  todopermutado[,j]=sample(submuestraribo[,j],71,replace=FALSE)
}



pithr=0.75

###############################################################################
#################################### graficos #################################

#grafico 1: stability seleccion en modelos graficos - dataset riboflavin
par(mfrow=c(2,5)) 
grafo(var(submuestraribo),0.27,TRUE)
grafo(var(submuestraribo),0.29,TRUE)
grafo(var(submuestraribo),0.3,TRUE)
grafo(var(submuestraribo),0.33,TRUE)
grafo(var(submuestraribo),0.35,TRUE)

stability.selection.grafico(submuestraribo,0.27,1000,floor(nrow(submuestraribo)/2),30,0.75,TRUE)
stability.selection.grafico(submuestraribo,0.29,1000,floor(nrow(submuestraribo)/2),30,0.75,TRUE)
stability.selection.grafico(submuestraribo,0.3,1000,floor(nrow(submuestraribo)/2),30,0.75,TRUE)
stability.selection.grafico(submuestraribo,0.33,1000,floor(nrow(submuestraribo)/2),30,0.75,TRUE)
stability.selection.grafico(submuestraribo,0.35,1000,floor(nrow(submuestraribo)/2),30,0.75,TRUE)


#grafico 2: stability seleccion en modelos graficos - todo permutado (grafo subyacente vacio)
par(mfrow=c(2,5)) 
grafo(var(todopermutado),0.08,TRUE)
grafo(var(todopermutado),0.07,TRUE)
grafo(var(todopermutado),0.06,TRUE)
grafo(var(todopermutado),0.05,TRUE)
grafo(var(todopermutado),0.04,TRUE)
stability.selection.grafico(todopermutado,0.08,100,floor(nrow(todopermutado)/2),30,0.75,TRUE)
stability.selection.grafico(todopermutado,0.07,100,floor(nrow(todopermutado)/2),30,0.75,TRUE)
stability.selection.grafico(todopermutado,0.06,100,floor(nrow(todopermutado)/2),30,0.75,TRUE)
stability.selection.grafico(todopermutado,0.05,100,floor(nrow(todopermutado)/2),30,0.75,TRUE)
stability.selection.grafico(todopermutado,0.04,100,floor(nrow(todopermutado)/2),30,0.75,TRUE)
