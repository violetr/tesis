######################### bibliotecas#########################

library(IsingFit)
library(here) # to manage finding files, avoiding harcoded personal paths
library(readr)
library(magrittr)
library(purrr)
library(tidyr)
library(dplyr)
library(stringr)
library(lubridate)
library(tidygraph)
library(ggraph)
library(extrafont)
loadfonts()

####################### cargo datos #########################

asuntos <- read_csv(here::here("datos/asuntos-diputados.csv"))
diputados <- read_csv(here::here("datos/diputados.csv"))
votos <- read_csv(here::here("datos/votaciones-diputados.csv"))
bloques_color <- read_csv(here::here("datos/bloques-diputados.csv"))

######################## ordeno #############################

asuntosselec <- asuntos %>%
  mutate(fecha = as.Date(fecha, "%m/%d/%Y")) %>%
  filter(fecha > ymd("2013-12-10"), fecha < ymd("2015-12-9")) %>%
  select(asuntoId)

votosselec <- votos %>%
  distinct(asuntoId, diputadoId, bloqueId, voto) %>%
  filter(
    asuntoId %in% asuntosselec$asuntoId,
    bloqueId %in% c(66, 67, 64, 109, 136, 179, 17, 18, 19, 172, 78),
    !(diputadoId %in% 1:10)
  )

bloques <- votosselec %>%
  select(diputadoId, bloqueId) %>%
  group_by(diputadoId) %>%
  top_n(1, bloqueId)%>%
  distinct(diputadoId, bloqueId)

pre_dataset <- votosselec %>%
  select(-bloqueId) %>%
  spread(key = diputadoId, value = voto) %>%
  select(-asuntoId)

# cambio labels
dataset <- pre_dataset
dataset[pre_dataset == 1 | pre_dataset == 2 | pre_dataset == 3] <- 0
dataset[pre_dataset == 0] <- 1

# Elimino los diputados que renunciaron,
# que faltaron mucho o sin varianza en el voto
# porque no los soporta IsingFit:
dataset <- dataset %>%
  select_if(~ all(!is.na(.))) %>% # los que no estuvieron el mandato completo?
  select_if(~ sum(.) < nrow(dataset) - 8) %>% # saco a los que votaron mucho +
  select_if(~ sum(.) > 8) # saco a los que faltaron mucho y o todo negativo

apellidos <- map(diputados[10:nrow(diputados), "nombre"], ~ str_split_fixed(.x, ",", n = 2))[[1]][, 1]

matching_tib <- tibble(apellidos = apellidos, ID = as.character(diputados[10:nrow(diputados), ]$diputadoID))

matching_dipu_id <- setNames(apellidos, as.character(diputados[10:nrow(diputados), ]$diputadoID))

diputados_colores <- left_join(bloques,
                               bloques_color,
                               by = c("bloqueId" = "id_bloque")) %>%
  ungroup() %>%
  mutate(diputadoId = as.character(diputadoId)) %>%
  left_join(matching_tib, by = c("diputadoId" = "ID")) %>%
  select(apellidos, color, bloqueId)%>%
  mutate(color= ifelse(is.na(color), '#85c1e9', color))

diputados_colores %>% filter(bloqueId==64)
unique(diputados_colores$apellidos)

matching_colores <- setNames(diputados_colores$color, diputados_colores$apellidos)
matching_bloques <- setNames(diputados_colores$bloqueId, diputados_colores$apellidos)
matching_bloques_colores <- setNames(unique(diputados_colores$color), unique(diputados_colores$bloqueId))

dipus_en_grafo <- matching_dipu_id[colnames(dataset)]

colnames(dataset) <- dipus_en_grafo

#############################################################

modelo <- IsingFit(dataset, gamma = 0.25)

#############################################################

grafo_crudo <- as_tbl_graph(modelo$weiadj, directed = FALSE)

grafo <- grafo_crudo %>%
  activate(nodes) %>%
  mutate(bloque = matching_bloques[name]) %>%
  activate(edges) %>%
  mutate(link = ifelse(weight > 0, 'positivo', 'negativo')) %>%
  mutate(weight = abs(weight))

pdf(here::here("figuras/pruebadipus36.pdf"), height = 15, width = 15)
ggraph(grafo) +
  geom_edge_link(aes(colour = color), alpha = 0.4) +
  scale_edge_colour_manual(values=c('#af270f', '#0faf62')) +
  geom_node_text(aes(label = name), repel = TRUE) +
  geom_node_point(aes(filter=bloque==67), size=3, color="#1f77b4" , alpha=0.5) +
  geom_node_point(aes(filter=bloque==172), size=3, color="#d62728" , alpha=0.5) +
  geom_node_point(aes(filter=bloque==78), size=3, color="#0b615e" , alpha=0.5) +
  geom_node_point(aes(filter=bloque==179), size=3, color="#e7ba52" , alpha=0.5) +
  geom_node_point(aes(filter=bloque==64), size=3, color="#85c1e9" , alpha=0.5) +
  theme_void()
dev.off()