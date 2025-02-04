
clades <- tibble(
  node =c(458, 453, 449, 448, 434, 431,
          331, 327, 342, 312, 315, 414,
          409, 384, 406, 378, 382, 380,
          405, 401, 397, 395, 388, 283,
          279, 294, 246, 247, 300, 306,
          308, 368, 371, 373, 365, 357,
          355, 350, 345, 343, 250, 252,
          253, 260, 262, 267, 268, 275),
  genus = c("Macaca", "Papio", "Cercocebus", "Mandrillus", "Cercopithecus", "Allochrocebus", "Alouatta", "Ateles", "Sapajus", "Cheracebus", "Plecturocebus", "Trachypithecus", "Semnopithecus", "Nomascus", "Pygathrix",  "Pongo", "Pan", "Gorilla", "Rhinopithecus", "Presbytis", "Piliocolobus", "Colobus", "Hylobates", "Eulemur", "Hapalemur", "Tarsius", "Galago", "Otolemur", "Pithecia", "Chiropotes", "Cacajao", "Callithrix", "Cebuella", "Mico", "Leontopithecus", "Saguinus", "Leontocebus", "Aotus", "Saimiri", "Cebus", "Perodicticus", "Loris", "Nycticebus", "Cheirogaleus", "Lepilemur", "Avahi", "Propithecus","Varecia"))

p1 <- p
for (i in 1:nrow(clades)){
  p1 <- collapse(p1, node = clades[i,]$node, mode = "max", clade_name = clades[i,]$genus, alpha = 1, color = "grey", fill = "light blue", size = 1)
  p1 <- p1 + geom_cladelab(node=clades[i,]$node, label=clades[i,]$genus, fontsize = 2, horizontal = FALSE, geom = "label", align = TRUE, extend = 1, offset = 1)
}

p1


collapsed_tree <-
  purrr::reduce(
    genNodes$node,
    \(x,y) collapse(x,y,mode = "max",fill="transparent",
                    color="black",size=0.1,),
    .init = p
  )

collapsed_tree +
  geom_cladelab(node=35, label="test label", angle=0,
                fontsize=8, offset=.5, vjust=.5)
