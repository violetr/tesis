#########################bibliotecas#########################

library(IsingFit)
library(here)
library(readr)
library(magrittr)
library(purrr)
library(lubridate)
library(tidyr)
library(dplyr)
library(stringr)
library(extrafont)
loadfonts()

####################### cargo datos #########################

asuntos = read_csv(here::here("datos/asuntos-diputados.csv"))
diputados = read_csv(here::here("datos/diputados.csv"))
votos = read_csv(here::here("datos/votaciones-diputados.csv"))
bloques_color = read_csv(here::here("datos/bloques-diputados.csv"))

######################## ordeno #############################

asuntosselec <- asuntos %>% 
  mutate(fecha = as.Date(fecha, "%m/%d/%Y"))%>%
  filter(fecha>ymd("2013-12-10"), fecha<ymd("2015-12-9"))%>%
  select(asuntoId)

votosselec <- votos %>%
  distinct(asuntoId, diputadoId, bloqueId, voto) %>%  
  filter(asuntoId %in% asuntosselec$asuntoId, 
         bloqueId %in% c(66, 67, 64, 109, 136, 179, 17, 18, 19, 172, 78),
         !(diputadoId %in% 1:10)) 

bloques <- votosselec %>%
  select(diputadoId,  bloqueId) %>%
  group_by(diputadoId) %>%
  top_n(1, bloqueId)

pre_dataset <- votosselec %>%
  select(-bloqueId) %>%
  spread(key = diputadoId, value = voto) %>%
  select(-asuntoId)

# cambio labels
dataset <- pre_dataset
dataset[pre_dataset==1 | pre_dataset==2 | pre_dataset==3] <- 0
dataset[pre_dataset==0] <- 1
# Elimino los diputados que renunciaron, 
# que faltaron mucho o sin varianza en el voto
# porque no los soporta IsingFit
dataset <- dataset %>%
  select_if(~all(!is.na(.))) %>%
  select_if(~sum(.) < nrow(dataset) - 8) %>%
  select_if(~sum(.) > 8)

apellidos <- map(diputados[10:nrow(diputados), "nombre"], ~ str_split_fixed(.x, ",", n= 2))[[1]][, 1]

matching_tib <- tibble(apellidos = apellidos, ID = as.character(diputados[10:nrow(diputados), ]$diputadoID))

matching_dipu_id <- setNames(apellidos, as.character(diputados[10:nrow(diputados), ]$diputadoID))

diputados_colores = left_join(bloques, 
                              bloques_color, 
                              by = c("bloqueId" = "id_bloque")) %>%
  ungroup()%>%
  mutate(diputadoId = as.character(diputadoId)) %>%
  left_join(matching_tib, by = c("diputadoId" = "ID")) %>%
  select(apellidos, color)

matching_colores = setNames(diputados_colores$color, diputados_colores$apellidos)

dipus_en_grafo <- matching_dipu_id[colnames(dataset)]

colnames(dataset) <- dipus_en_grafo

#############################################################

modelo = IsingFit(dataset, gamma = 0.25)

ady = modelo$weiadj !=0 
for (i in 1:dim(ady)[1]) {
  ady[i, i] = FALSE
}

negativos<-1*(modelo$weiadj<0)
pesos<- abs(modelo$weiadj)
length(negativos)

# borro los aislados
sinisolatedd = ady[apply(ady, 2, sum, na.rm=TRUE) != 0, apply(ady, 2, sum, na.rm=TRUE) != 0]
bloquessiniso = bloques[apply(ady, 2, sum) != 0]
negativossiniso <- negativos[apply(ady, 2, sum) != 0,apply(ady, 2, sum) != 0]
pesossiniso <- pesos[apply(ady, 2, sum) != 0,apply(ady, 2, sum) != 0]
#####################colores#################################

#coloresh <- c("firebrick2", "dodgerblue2", "gold", "firebrick1", "deepskyblue1", "gold")
coloresh <- c("paleturquoise1","gold","indianred1","lightgrey", "skyblue")

coloresh2 <- c("red","blue", "orange", "green", "red", "yellow", "paleturquoise1", "blue", "black", "red","black","black", "red","red", "blue", "black", "black", "black", "blue", "black", "black", "black","black","red")
bloques.df[bloques.df$bloqueId == 86, 2]
names(coloresh2) = unique(bloques)

#indianred1
# 172 ucr firebrick
# 67  FPV dodgerblue2
# 136 pro gold
# 18 coalicion civica - ARI firebrick1
# 109 nuevo encuentro dodgerblue2
# 64 frente nuevo encuentro
# 179 union pro gold

##############################################################

graph <- graph_from_adjacency_matrix(sinisolatedd)

# Not specifying the layout - defaults to "auto"


netd <- network::network(sinisolatedd, directed = FALSE)

netd %v% "bloque" <- bloquessiniso

netd %e% "negativo" <- negativossiniso

netd %e% "weights" <- pesossiniso

pdf("newwwwww3.pdf",height = 6, width = 18)
ggnet2(netd, 
       size = 7, 
       label = TRUE, 
       alpha = 1, 
       label.size = 2, 
       color = coloresh2[netd %v% "bloque"], 
       edge.color = ifelse(netd %e% "negativo"==1,"firebrick1","forestgreen"))
dev.off()

netd <- network::network(ady, directed = FALSE)

netd %v% "bloque" <- bloques

netd %e% "negativo" <- negativos



pdf("holuuuuu68.pdf",height = 6, width = 18)
ggnet2(netd, 
       size = 7, 
       label = TRUE, 
       alpha = 1, 
       label.size = 2, 
       color = coloresh[netd %v% "bloque"], 
       edge.color = ifelse(netd %e% "negativo"==1,"firebrick1","forestgreen"))
dev.off()

################################################################################

holu = igraph::graph_from_adjacency_matrix(sinisolatedd, weighted=TRUE)
holuu = as_tbl_graph(holu, directed = FALSE)

holii = as_tbl_graph(modelo$weiadj, directed = FALSE)

holuu = holii %>% 
  activate(nodes) %>%
  mutate(color_bloque = matching_colores[name]) %>%
  activate(edges) %>% 
  mutate(color = weight<0) %>%
  mutate(weight = abs(weight))

holuu

pdf(here::here("figuras/pruebadipus14.pdf"), height = 15, width = 15)
ggraph(holuu, layout = 'kk') + 
  geom_edge_link()+
  geom_node_point(aes(colour=color_bloque))+ theme_void()+
  geom_node_text(aes(label = name), repel = TRUE)
dev.off()

pdf(here::here("figuras/pruebadipus1.pdf"), height = 12, width = 12)
ggraph(holuu,layout="mds") + 
  geom_edge_link(aes(width = abs(weight), colour = color), alpha = 0.8) + 
  geom_node_point(aes(colour=color.bloque))+ theme_graph()+
  geom_node_text(aes(label = name), repel = TRUE)
dev.off()


