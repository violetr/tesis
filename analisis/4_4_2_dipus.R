#########################bibliotecas#########################


if (! require(IsingFit)) {
  install.packages(IsingFit)
  library(IsingFit)
}

source('./metodos/seleccion_modelo_normal.R')
source('./auxiliares.R')

#############################################################

#auxiliares especificos para este dataset

dameIndiceAsunto <- function(asuntoId, asuntosselec) {
  return(which(asuntosselec == asuntoId))
}

dameListaIndicesAsuntos <- function(asuntosId, asuntosselec) {
  lista = vector()
  for (i in 1:length(asuntosId)) {
    lista[i] = dameIndiceAsunto(asuntosId[i], asuntosselec)
  }
  return(lista)
}

dameApellidoDiputado <- function(DiputadoId, datadiputados) {
  return(unlist(strsplit(as.character(diputados[DiputadoId,2]), split = ", "))[1])  
}

dameListaApellidos <- function(conjDiputadoId, datadiputados) {
  lista = vector()
  for(i in 1:length(conjDiputadoId)){
    lista[i] = dameApellidoDiputado(conjDiputadoId[i], datadiputados)
  }
  return(lista)
}

dameBloque <- function(DiputadoId, datavotos){
  return(datavotos[datavotos$diputadoId == DiputadoId, 3][1])
}

dameListaBloques <- function(conjDiputadoId, datavotos){
  lista = vector()
  for (i in 1:length(conjDiputadoId)) {
    lista[i] = dameBloque(conjDiputadoId[i], datavotos)
  }
  return(lista)
}

contarFaltas <- function(columna) {
  return(sum(is.na(columna)))
}


############################datos##############################

#as.Date(asuntos$fecha,format='%m/%d/%Y')<as.Date('05/12/2002',format='%m/%d/%Y')
#  as.Date("2007-06-22")


votos <- read.csv('./datos/votaciones-diputados.csv')
diputados <- read.csv('./datos/diputados.csv')
asuntos <- read.csv('./datos/asuntos-diputados.csv')
asuntos$fecha

#asuntosselec = asuntos[asuntos$ano >= 2013 & asuntos$ano <= 2015, 1]
esta.en.tiempo = as.Date(asuntos$fecha,format='%m/%d/%Y')>as.Date('10/12/2013',format='%m/%d/%Y') & as.Date(asuntos$fecha,format='%m/%d/%Y')<as.Date('9/12/2015',format='%m/%d/%Y')
asuntosselec = asuntos[esta.en.tiempo, 1]




length(asuntosselec)
# me quedo con los votos corresp a los asuntos seleccionados cuyos bloques pertenecen
# a los mas importantes y que no sean el diputado 5 que no esta identificado
votosselec = votos[votos$asuntoId  %in% asuntos[esta.en.tiempo, 1] & votos$bloqueId %in% c(66, 67, 64, 109, 136, 179, 17, 18, 19, 172) & votos$diputadoId != 5, ]

diputadosselec = unique(votosselec$diputadoId)
#borro a un par que renunciaron y a los que los reemplazaron
diputadosselec <- setdiff(diputadosselec, c(296,912,1037,753))
#dameListaApellidos(diputadosselec,diputados)
#verificar carac
#length(diputadosselec)
#unique(votosselec[votosselec$diputadoId %in% diputadosselec, "bloqueId"])
#unique(dameListaBloques(diputadosselec,votosselec))
#unique(votosselec$bloqueId) # hay uno que pertenecio a dos bloques

nasuntos = length(asuntosselec)
ndiputados = length(diputadosselec)

mis.datitos = as.data.frame(matrix(rep(0, nasuntos * ndiputados), nrow = nasuntos))

# negativo, abstencion y ausente seran computados como 0, positivo como 1.
# primero miro las abstenciones para eliminar a los que faltaron mucho
for (j in 1:length(diputadosselec)) {
  #busco los votos afirmativos (factor 0) y los pongo como 1 en la matriz de datos el valor 0 corresponde
  # a voto neg abstencion y ausencia
  votosafirmativos = votosselec[votosselec$diputadoId == diputadosselec[j] & votosselec$voto == 0, ]$asuntoId
  ausencias = votosselec[votosselec$diputadoId == diputadosselec[j] & votosselec$voto == 3, ]$asuntoId
  nueva.columna = rep(0, nasuntos)
  indicesAfirmativos = dameListaIndicesAsuntos(votosafirmativos, asuntosselec)
  indicesAusencias = dameListaIndicesAsuntos(ausencias, asuntosselec)
  nueva.columna[indicesAfirmativos] = 1
  nueva.columna[indicesAusencias] = NA
  mis.datitos[ , j] = nueva.columna
}

nrosvar = obtenerindices(apply(mis.datitos, 2, contarFaltas) < nrow(mis.datitos) /4 & apply(mis.datitos, 2, sum,na.rm=TRUE) > 8 & (nrow(mis.datitos) - apply(mis.datitos, 2, sum,na.rm=TRUE)) > 8)


nuevadata = mis.datitos[ , nrosvar]
nuevadata[is.na(nuevadata)] = 0

dim(nuevadata)
nombrescolumnas = dameListaApellidos(diputadosselec[nrosvar], diputados)
colnames(nuevadata) = nombrescolumnas

bloques = as.character(dameListaBloques(diputadosselec[nrosvar], votosselec))
#############################################################

modelo = IsingFit(nuevadata, gamma = 0.18)

names(modelo)

plot(modelo)

ady = modelo$weiadj !=0 
for (i in 1:dim(ady)[1]) {
  ady[i, i] = FALSE
}

negativos<-1*(modelo$weiadj<0)


# borro los aislados
sinisolatedd = ady[apply(ady, 2, sum) != 0, apply(ady, 2, sum) != 0]
dim(sinisolatedd)
bloquessiniso = bloques[apply(ady, 2, sum) != 0]
negativossiniso <- negativos[apply(ady, 2, sum) != 0,apply(ady, 2, sum) != 0]

#####################colores#################################

#coloresh <- c("firebrick2", "dodgerblue2", "gold", "firebrick1", "deepskyblue1", "gold")
coloresh <- c("paleturquoise1","indianred1", "skyblue", "gold","gold")

names(coloresh) = unique(bloques)

#indianred1
# 172 ucr firebrick
# 67  FPV dodgerblue2
# 136 pro gold
# 18 coalicion civica - ARI firebrick1
# 109 nuevo encuentro dodgerblue2
# 179 union pro gold

##############################################################

netd <- network::network(sinisolatedd, directed = FALSE)

netd %v% "bloque" <- bloquessiniso

netd %e% "negativo" <- negativossiniso

ggnet2(netd, size = 7, label = TRUE, alpha = 1, label.size = 2, color = coloresh[netd %v% "bloque"], edge.color = ifelse(netd %e% "negativo"==1,"firebrick1","forestgreen"))

