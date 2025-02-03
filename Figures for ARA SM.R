library(tidyverse)
library(treedata.table)
library(ape)
library(ggtree)
library(ggrepel)
library(castor)
library(phytools)
library(phangorn)
library(here)
library(readxl)

rm(list = ls())

colors <- c("blue", "green", "orange", "red","maroon","skyblue")
state_names <- c("G", "P", "S")

# Plot SM Figure 2 ----
tree_file <- "Kuderna data_s4_fossil_calibrated_time_tree.nex.tree"
tree <- read.tree(tree_file)
outgroup <- c("Tupaia_belangeri", "Galeopterus_variegatus", "Mus_musculus", "Oryctolagus_cuniculus")
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)

# quick plot
plot.phylo(tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
# note: tips are not aligned

d <- read_csv("Kuderna et al taxa.csv", col_names = TRUE) # adds in other taxonomic levels for taxa in Kuderna et al data set
d <- d |>
  rowwise() |>
  mutate(Genus = strsplit(label, "_")[[1]][1])

# make treedata.table object
t <- as.treedata.table(tree = tree, data = as.data.frame(d))

# t$dat holds the data...
t$dat <- t$dat |>
  rowwise() |>
  mutate(Genus = strsplit(tip.label, "_")[[1]][1])

# t$phy holds the phylogeny
tree <- t$phy

# force the tree to be ultrametric for visualization
um_tree <- force.ultrametric(tree)
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
# now tips are aligned

# plot with ggplot
p <- ggtree(um_tree, size = 0.3, layout = "fan") %<+% d + geom_tiplab(aes(label=label), size=2)

# code below is to get nodes for MRCA of each Superfamily and make dataframe for clade nodes
clades <- tibble(
  clade = c(unique(d$Superfamily)),
  node = NA
) |>
  filter(!is.na(clade) & clade != "Primates")

# find the MRCA for each clade
for (i in 1:length(clades$clade)) {
  clades$node[i] <- MRCA(
    tree,
    d |> filter(Superfamily == clades$clade[i]) |> pull(label)
  )
}

# code below is to get nodes for MRCA of each Genus and make dataframe for clade nodes
genus_clades <- tibble(
  clade = c(unique(d$Genus)),
  node = NA
) |>
  filter(!is.na(clade))

# find the MRCA for each clade
for (i in 1:length(genus_clades$clade)) {
  genus_clades$node[i] <- MRCA(
    tree,
    d |> filter(Genus == genus_clades$clade[i]) |> pull(label)
  )
}

# final plot SM Figure 2

p <- ggtree(um_tree, size = 0.3, layout = "fan") %<+% d +
  # geom_tiplab(aes(label=label), size=2) +
  geom_highlight(data=clades,
               aes(node=node,
                   fill=clade),
               alpha=0.5,
               align="left") + #,
               # show.legend = FALSE) +
  scale_fill_manual(
    values=colors,
    breaks = c("Lemuroidea", "Lorisoidea", "Tarsioidea", "Ceboidea", "Hominoidea", "Cercopithecoidea")) +
  geom_cladelab(data=genus_clades,
                mapping=aes(node=node, label=clade),
                fontsize=3.5,
                align="TRUE",
                angle="auto",
                offset=1,
                offset.text=1) +
  labs(fill = "Primate Superfamily") +
  geom_tree(linewidth=0.3) +
  geom_tippoint(size = 0.2) +
  theme(legend.position="left",
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=12)) #change legend text font size)

p # save this as a high resolution PNG file

rm(list = ls())

# Process Olivier et al dataset to generate base data ----
## Load character data from Olivier et al, with additional taxonomic levels added
base_data <- read_csv("Olivier et al 2024 data/dataG.csv", col_names = TRUE)
base_data <- base_data |>
  rowwise() |>
  mutate(`Genus` = str_split(`Genus_species`, "_")[[1]][1]) |>
  ungroup() |>
  rename(`Species` = `Genus_species`) |>
  mutate(`main_SO` = `Main1`) |>
  filter(`sp.pop` == "sp") |>
  mutate(
    character_data = case_when(
      main_SO == "solitary" ~ "S",
      main_SO == "MF" ~ "P",
      main_SO == "MFF" ~ "G",
      main_SO == "FFMM" ~ "G",
      main_SO == "FMM" ~ "G",
      main_SO == "MF_MFF" ~ "G"
      )
    ) |>
  select(Order, Superfamily, Genus, Species, character_data)

## Save this base data file
write_csv(base_data, "base_data.csv")

rm(list = ls())

# Data wrangling to add to base dataset ----

## Load in updated Olivier et al data... we keep these base data for all primates other than nocturnal strepsirrhines and tarsiers
base_data <- read_csv("base_data.csv", col_names = TRUE)

## Columns to keep
keep_columns <- c("Order", "Superfamily", "Family", "Species", "Kappeler & Pozzi 3 state","Lukas & Clutton-Brock 3 state", "Olivier et al 3 state", "Shultz et al 3 state", "Müller & Thalman 3 state", "Müller & Thalman 3 state - dispersed", "Additional 3 state", "Additional 3 state - dispersed", "Conflicting Character State Assignment?")

## Load in comparison dataset from Excel... we join this to base data from Olivier et al for all primates other than nocturnal strepsirrhines and tarsiers
comparison_data <- read_excel("comparison table.xlsx", sheet = 1, col_names = TRUE) |>
  select(all_of(keep_columns)) |>
  mutate(`Species` = str_replace(`Species`, " ", "_"))

combined_data <- full_join(base_data, comparison_data, by = c("Species" = "Species")) |>
  mutate(Order.y = if_else(is.na(Order.y), Order.x, Order.y)) |>
  mutate(Order.x = if_else(is.na(Order.x), Order.y, Order.x)) |>
  mutate(Superfamily.y = if_else(is.na(Superfamily.y), Superfamily.x, Superfamily.y)) |>
  mutate(Superfamily.x = if_else(is.na(Superfamily.x), Superfamily.y, Superfamily.x)) |>
  select(-c("Order.y", "Superfamily.y", "Family")) |> rename(`Order` = `Order.x`, `Superfamily` = `Superfamily.x`) |>
  distinct() |>
  mutate(keep_index = !is.na(`character_data`) & !(Superfamily %in% c("Lemuroidea", "Lorisoidea", "Tarsioidea"))) |>
  rowwise() |>
  mutate(keep_index = if_else(keep_index == FALSE & Genus %in% c("Lemur", "Eulemur", "Hapalemur", "Prolemur", "Varecia", "Propithecus", "Indri"), TRUE, keep_index)) |>
  mutate(`Olivier et al 3 state` = if_else(keep_index == TRUE & is.na(`Olivier et al 3 state`), `character_data`, `Olivier et al 3 state`)) |>
  mutate(`Kappeler & Pozzi 3 state` = if_else(keep_index == TRUE & is.na(`Kappeler & Pozzi 3 state`), `character_data`, `Kappeler & Pozzi 3 state`)) |>
  mutate(`Lukas & Clutton-Brock 3 state` = if_else(keep_index == TRUE & is.na(`Lukas & Clutton-Brock 3 state`), `character_data`, `Lukas & Clutton-Brock 3 state`)) |>
  mutate(`Shultz et al 3 state` = if_else(keep_index == TRUE & is.na(`Shultz et al 3 state`), `character_data`, `Shultz et al 3 state`)) |>
  mutate(`Müller & Thalman 3 state` = if_else(keep_index == TRUE & is.na(`Müller & Thalman 3 state`), `character_data`, `Müller & Thalman 3 state`)) |>
  mutate(`Müller & Thalman 3 state - dispersed` = if_else(keep_index == TRUE & is.na(`Müller & Thalman 3 state - dispersed`), `character_data`, `Müller & Thalman 3 state - dispersed`)) |>
  select(-c("character_data", "keep_index"))

## Save this data file
write_csv(combined_data, "combined_data.csv")

rm(list = ls())

combined_data <- read_csv("combined_data.csv", col_names = TRUE)

## Now create various datasets based on combined base_data + nocturnal strepsirrhine, tarsier, and outgroup character states
d <- list() # list to hold datasets
s <- list() # list to hold study names

## Olivier et al
d[[1]] <- combined_data |>
  filter(!is.na(`Olivier et al 3 state`)) |>
  mutate(character_data = `Olivier et al 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[1]] <- "Olivier et al"

## Olivier et al - adding in my data for outgroup taxa
d[[2]] <- combined_data |>
  filter(!is.na(`Olivier et al 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Olivier et al 3 state`), `Olivier et al 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[2]] <- "Olivier et al additional"

## Kappeler & Pozzi
d[[3]] <- combined_data |>
  filter(!is.na(`Kappeler & Pozzi 3 state`)) |>
  mutate(character_data = `Kappeler & Pozzi 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[3]] <- "Kappeler & Pozzi"

## Kappeler & Pozzi - adding in my data for outgroup taxa
d[[4]] <- combined_data |>
  filter(!is.na(`Kappeler & Pozzi 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Kappeler & Pozzi 3 state`), `Kappeler & Pozzi 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[4]] <- "Kappeler & Pozzi additional"

## Lukas & Clutton-Brock
d[[5]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`)) |>
  mutate(character_data = `Lukas & Clutton-Brock 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[5]] <- "Lukas & Clutton-Brock"

## Lukas & Clutton-Brock - adding in my data for outgroup taxa they didn't have and replacing theirs that is incorrect
d[[6]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Additional 3 state`), `Additional 3 state`, `Lukas & Clutton-Brock 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[6]] <- "Lukas & Clutton-Brock additional"

## Lukas & Clutton-Brock - adding in my data only for outgroup taxa they didn't have
d[[11]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Lukas & Clutton-Brock 3 state`), `Lukas & Clutton-Brock 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[11]] <- "Lukas & Clutton-Brock additional 2"

## Müller & Thalman
d[[7]] <- combined_data |>
  filter(!is.na(`Müller & Thalman 3 state`)) |>
  mutate(character_data = `Müller & Thalman 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[7]] <- "Müller & Thalman"

## Müller & Thalman - adding in my data for outgroup taxa
d[[8]] <- combined_data |>
  filter(!is.na(`Müller & Thalman 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Müller & Thalman 3 state`), `Müller & Thalman 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[8]] <- "Müller & Thalman additional"

## Müller & Thalman - dispersed
d[[9]] <- combined_data |>
  filter(!is.na(`Müller & Thalman 3 state - dispersed`)) |>
  mutate(character_data = `Müller & Thalman 3 state - dispersed`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "D-P, D-G", "D-P", `character_data`))
s[[9]] <- "Müller & Thalman - dispersed"

## Müller & Thalman - dispersed - adding in my data for outgroup taxa
d[[10]] <- combined_data |>
  filter(!is.na(`Müller & Thalman 3 state - dispersed`) | !is.na(`Additional 3 state - dispersed`)) |>
  mutate(character_data = if_else(!is.na(`Müller & Thalman 3 state - dispersed`), `Müller & Thalman 3 state - dispersed`, `Additional 3 state - dispersed`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "D-P, D-G", "D-P", `character_data`))
s[[10]] <- "Müller & Thalman additional - dispersed"

rm(list = setdiff(ls(), c("d", "s")))

# Analyses on Kuderna et al tree ----

## Get full tree and plot ----
tree_file <- "Kuderna data_s4_fossil_calibrated_time_tree.nex.tree"
tree <- read.tree(tree_file)

# Replace Cephalopachus and Carlito with Tarsis
tree$tip.label <- gsub("Cephalopachus", "Tarsius", tree$tip.label)
tree$tip.label <- gsub("Carlito", "Tarsius", tree$tip.label)

outgroup <- c("Mus_musculus", "Oryctolagus_cuniculus", "Tupaia_belangeri", "Galeopterus_variegatus")
is.rooted.phylo(tree)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
Kuderna_et_al_tree <- tree

um_tree <- force.ultrametric(tree, method = "extend")

# plot full tree
p <- ggtree(um_tree,
            layout = "fan",
            size = 0.3,
            branch.length = "branch.length"
) +
  geom_tiplab(size = 2)
p

## Repeat for each dataset ----
Kuderna_et_al_tree_res <- tibble(dataset = character(), ntaxa = numeric(), ntips = numeric(), phylo = character(), method1 = character(), model1 = character(), DG1 = numeric(), DP1 = numeric(), G1 = numeric(), P1 = numeric(), S1 = numeric(), sum1 = numeric(), method2 = character(), model2 = character(), DG2 = numeric(), DP2 = numeric(), G2 = numeric(), P2 = numeric(), S2 = numeric(), sum2 = numeric())

for (i in 1:length(d)){
  if (i == 9 | i == 10){
    colors <- c("skyblue", "red", "blue", "green", "orange")
    state_names <- c("D-G", "D-P", "G", "P", "S")
  } else {
    colors <- c("blue", "green", "orange", "red","maroon","skyblue")
    state_names <- c("G", "P", "S")
  }
  ## Create treedata.table to merge tree and data and drop tips that are missing in dataset and data that are missing in tree
  tree <- Kuderna_et_al_tree
  t <- as.treedata.table(tree = tree, data = as.data.frame(d[[i]]), name_column = "Species")

  phy <- t$phy # get the phylogeny
  is.rooted.phylo(phy) # check if it's rooted
  is.binary(phy) # check if it's binary
  outgroup <- t$dat |> # set outgroup to be non primates if they are present in the dataset
    filter(`Order` %in% c("Lagomorpha", "Scandentia", "Rodentia", "Dermoptera")) |>
    pull(`tip.label`)
  if (length(outgroup) == 0) { # otherwise set outgroups to be strepsirrhines
    outgroup <- t$dat |>
      filter(`Superfamily` %in% c("Lemuroidea", "Lorisoidea")) |>
      pull(`tip.label`)
  }

  t$phy <- root(t$phy, outgroup = outgroup, resolve.root = TRUE) # reroot the tree
  phy <- t$phy
  is.rooted.phylo(phy) # check if it's rooted
  is.binary(phy) # check if it's binary

  ### Run ASR -----
  #### Set Up Dataset ----
  character_data <- t$dat$`character_data`
  names(character_data) <- t$dat$tip.label
  character_data <- as.factor(character_data) # convert states to a factor (required for discrete traits in phytools)
  tree <- t$phy

  #### Run ML Ancestral State Reconstruction with ER ----
  ace_results <- ace(character_data, tree, type = "discrete", method = "ML", model = "ER") # default = marginal = FALSE, returns empirical Bayesian posterior probabilities

  #### Run SCM Ancestral State Reconstruction with ER ----
  n_simulations <- 100  # Number of stochastic maps to generate
  scm_results <- make.simmap(tree, character_data, model = "ER", nsim = n_simulations, pi = "fitzjohn")

  # Summarize the stochastic maps
  summary_maps <- summary(scm_results)

  # Extract posterior probabilities and pies for internal nodes
  posterior_probs <- summary_maps$ace
  sim_node_pies <- posterior_probs[1:tree$Nnode, ]

  # Tip states
  sim_tip_pies <- summary_maps$tips

  #### Plotting ML ----
  um_tree <- force.ultrametric(tree, method = "extend")
  um_tree$root.edge <- 2

  plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)

  # Add pie charts for ancestral states at internal nodes
  node_pies <- as.matrix(ace_results$lik.anc)
  nodelabels(
    pie = node_pies,
    piecol = colors,
    cex = 0.2
  )

  # Add pie charts for tips
  tip_pies <- as.factor(character_data)
  tip_pies <- to.matrix(tip_pies, levels(tip_pies))
  tiplabels(
    pie = tip_pies,
    piecol = colors,
    cex = 0.2
  )

  legend("topleft",
         legend = state_names,
         fill = colors,
         title = "Character States")

  mrca <- MRCA(
    t$phy,
    t$dat |> filter(Order == "Primates") |> pull(tip.label))

  root_pie <- subset(node_pies, rownames(node_pies) %in% mrca)
  sim_root_pie <- subset(sim_node_pies, rownames(sim_node_pies) %in% mrca)

  if (i == 9 | i == 10){
    r <- tibble(dataset = s[[i]], phylo = "Kuderna et al", ntaxa = nrow(t$dat), ntips = length(t$phy$tip.label), method1 = "ML", model1 = "ER", DG1 = root_pie[,"D-G"], DP1 = root_pie[,"D-P"], G1 = root_pie[,"G"], P1 = root_pie[,"P"], S1 = root_pie[,"S"], sum1 = DG1 + DP1 + G1 + P1 + S1, method2 = "SCM", model2 = "ER", DG2 = sim_root_pie[,"D-G"], DP2 = sim_root_pie[,"D-P"], G2 = sim_root_pie[, "G"], P2 = sim_root_pie[ ,"P"], S2 = sim_root_pie[ ,"S"], sum2 = DG2 + DP2 + G2 + P2 + S2)
  } else {
    r <- tibble(dataset = s[[i]], phylo = "Kuderna et al", ntaxa = nrow(t$dat), ntips = length(t$phy$tip.label), method1 = "ML", model1 = "ER", DG1 = 0, DP1 = 0, G1 = root_pie[,"G"], P1 = root_pie[,"P"], S1 = root_pie[,"S"], sum1 = DG1 + DP1 + G1 + P1 + S1, method2 = "SCM", model2 = "ER", DG2 = 0, DP2 = 0, G2 = sim_root_pie[, "G"], P2 = sim_root_pie[ ,"P"], S2 = sim_root_pie[ ,"S"], sum2 = DG2 + DP2 + G2 + P2 + S2)
  }

  Kuderna_et_al_tree_res <- bind_rows(Kuderna_et_al_tree_res, r)
}

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree")))

# Analyses on Olivier et al tree ----

## Get full tree and plot ----
# start with multiple to capture uncertainty
tree_file <- "Olivier et al 2024 data/vert phylo.nex"
trees = read.nexus(tree_file)
base_data <- read_csv("base_data.csv", col_names = TRUE)

# create consensus phylogeny for robustness checks
tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")
is.rooted.phylo(tree)
outgroup <- base_data |> filter(Superfamily %in% c("Lemuroidea", "Lorisoidea")) |> pull(Species)

tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
tree$node.label <- NULL
Olivier_et_al_tree <- tree

um_tree <- force.ultrametric(tree, method = "extend")

# plot full tree
p <- ggtree(um_tree,
            layout = "circular",
            size = 0.3,
            branch.length = "branch.length"
) +
  geom_tiplab(size = 2)
p

## Repeat for each dataset ----
Olivier_et_al_tree_res <- tibble(dataset = character(), ntaxa = numeric(), ntips = numeric(), phylo = character(), method1 = character(), model1 = character(), DG1 = numeric(), DP1 = numeric(), G1 = numeric(), P1 = numeric(), S1 = numeric(), sum1 = numeric(), method2 = character(), model2 = character(), DG2 = numeric(), DP2 = numeric(), G2 = numeric(), P2 = numeric(), S2 = numeric(), sum2 = numeric())

for (i in 1:length(d)){
  if (i == 9 | i == 10){
    colors <- c("skyblue", "red", "blue", "green", "orange")
    state_names <- c("D-G", "D-P", "G", "P", "S")
  } else {
    colors <- c("blue", "green", "orange", "red","maroon","skyblue")
    state_names <- c("G", "P", "S")
  }
  ## Create treedata.table to merge tree and data and drop tips that are missing in dataset and data that are missing in tree
  tree <- Olivier_et_al_tree
  t <- as.treedata.table(tree = tree, data = as.data.frame(d[[i]]), name_column = "Species")

  phy <- t$phy # get the phylogeny
  is.rooted.phylo(phy) # check if it's rooted
  is.binary(phy) # check if it's binary
  outgroup <- t$dat |> # set outgroup to be non primates if they are present in the dataset
    filter(`Order` %in% c("Lagomorpha", "Scandentia", "Rodentia", "Dermoptera")) |>
    pull(`tip.label`)
  if (length(outgroup) == 0) { # otherwise set outgroups to be strepsirrhines
    outgroup <- t$dat |>
      filter(`Superfamily` %in% c("Lemuroidea", "Lorisoidea")) |>
      pull(`tip.label`)
  }

  t$phy <- root(t$phy, outgroup = outgroup, resolve.root = TRUE) # reroot the tree
  phy <- t$phy
  is.rooted.phylo(phy) # check if it's rooted
  is.binary(phy) # check if it's binary

  ### Run ASR -----
  #### Set Up Dataset ----
  character_data <- t$dat$`character_data`
  names(character_data) <- t$dat$tip.label
  character_data <- as.factor(character_data) # convert states to a factor (required for discrete traits in phytools)
  tree <- t$phy
  tree$node.label <- 1:Nnode(tree) + length(tree$tip.label)

  #### Run ML Ancestral State Reconstruction with ER ----
  ace_results <- ace(character_data, tree, type = "discrete", method = "ML", model = "ER") # default = marginal = FALSE, returns empirical Bayesian posterior probabilities

  #### Run SCM Ancestral State Reconstruction with ER ----
  n_simulations <- 100  # Number of stochastic maps to generate
  scm_results <- make.simmap(tree, character_data, model = "ER", nsim = n_simulations, pi = "fitzjohn")

  # Summarize the stochastic maps
  summary_maps <- summary(scm_results)

  # Extract posterior probabilities and pies for internal nodes
  posterior_probs <- summary_maps$ace
  sim_node_pies <- posterior_probs[1:tree$Nnode, ]

  # Tip states
  sim_tip_pies <- summary_maps$tips

  #### Plotting ML ----
  um_tree <- force.ultrametric(tree, method = "extend")
  um_tree$root.edge <- 2

  plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)

  # Add pie charts for ancestral states at internal nodes
  node_pies <- as.matrix(ace_results$lik.anc)
  nodelabels(
    pie = node_pies,
    piecol = colors,
    cex = 0.2
  )

  # Add pie charts for tips
  tip_pies <- as.factor(character_data)
  tip_pies <- to.matrix(tip_pies, levels(tip_pies))
  tiplabels(
    pie = tip_pies,
    piecol = colors,
    cex = 0.2
  )

  legend("topleft",
         legend = state_names,
         fill = colors,
         title = "Character States")

  mrca <- MRCA(
    t$phy,
    t$dat |> filter(Order == "Primates") |> pull(tip.label))

  root_pie <- subset(node_pies, rownames(node_pies) %in% mrca)
  sim_root_pie <- subset(sim_node_pies, rownames(sim_node_pies) %in% mrca)

  if (i == 9 | i == 10){
    r <- tibble(dataset = s[[i]], phylo = "Olivier et al", ntaxa = nrow(t$dat), ntips = length(t$phy$tip.label), method1 = "ML", model1 = "ER", DG1 = root_pie[ ,"D-G"], DP1 = root_pie[ ,"D-P"], G1 = root_pie[ ,"G"], P1 = root_pie[ , "P"], S1 = root_pie[, "S"], sum1 = DG1 + DP1 + G1 + P1 + S1, method2 = "SCM", model2 = "ER", DG2 = sim_root_pie[ ,"D-G"], DP2 = sim_root_pie[ ,"D-P"], G2 = sim_root_pie[ , "G"], P2 = sim_root_pie[ , "P"], S2 = sim_root_pie[ ,"S"], sum2 = DG2 + DP2 + G2 + P2 + S2)
  } else {
    r <- tibble(dataset = s[[i]], phylo = "Olivier et al", ntaxa = nrow(t$dat), ntips = length(t$phy$tip.label), method1 = "ML", model1 = "ER", DG1 = 0, DP1 = 0, G1 = root_pie[ ,"G"], P1 = root_pie[ , "P"], S1 = root_pie[, "S"], sum1 = DG1 + DP1 + G1 + P1 + S1, method2 = "SCM", model2 = "ER", DG2 = 0 , DP2 = 0, G2 = sim_root_pie[ , "G"], P2 = sim_root_pie[ , "P"], S2 = sim_root_pie[ ,"S"], sum2 = DG2 + DP2 + G2 + P2 + S2)
  }

  Olivier_et_al_tree_res <- bind_rows(Olivier_et_al_tree_res, r)
}

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree")))

# Analyses of effects of taxon sampling ----
# Here we can filter for taxa of interest or take random samples

root_pies <- tibble(dataset = numeric(), rep = numeric(), primateMRCA = numeric(), DG = numeric(), DP = numeric(), G = numeric(), P = numeric(), S = numeric())

colors <- c("skyblue", "red", "blue", "green", "orange")
state_names <- c("D-G", "D-P", "G", "P", "S")

# Generate 100 samples with 1 species per genus for each dataset and do ASR

for (i in 1:length(d)){
  for (j in 1:100) {
    for (k in c(Kuderna_et_al_tree, Olivier_et_al_tree)){
    tree <- k
    # pick one species from each Genus where the Species is represented in the given phylogeny
    data_subset <- d[[i]] |>
      group_by(Genus) |>
      filter(`Species` %in% tree$tip.label) |>
      sample_n(size = 1, replace = FALSE) |> # sample 1 species per genus
      pull(`Species`)

    t <- as.treedata.table(tree = tree, data = as.data.frame(d[[i]]), name_column = "Species")

    t <- t[tip.label %in% data_subset,]

    phy <- t$phy # get the phylogeny
    is.rooted.phylo(phy) # check if it's rooted
    is.binary(phy) # check if it's binary
    outgroup <- t$dat |> # set outgroup to be non primates if they are present in the dataset
      filter(`Order` %in% c("Lagomorpha", "Scandentia", "Rodentia", "Dermoptera")) |>
      pull(`tip.label`)
    if (length(outgroup) == 0) { # otherwise set outgroups to be strepsirrhines
      outgroup <- t$dat |>
        filter(`Superfamily` %in% c("Lemuroidea", "Lorisoidea")) |>
        pull(`tip.label`)
    }

    t$phy <- root(t$phy, outgroup = outgroup, resolve.root = TRUE) # reroot the tree
    phy <- t$phy
    is.rooted.phylo(phy) # check if it's rooted
    is.binary(phy) # check if it's binary

    ### Run ASR -----
    #### Set Up Dataset ----
    character_data <- t$dat$`character_data`
    names(character_data) <- t$dat$tip.label
    character_data <- as.factor(character_data) # convert states to a factor (required for discrete traits in phytools)
    tree <- t$phy

    #### Run ML Ancestral State Reconstruction----
    ace_results <- ace(character_data, tree, type = "discrete", method = "ML", model = "ER")

    #### Plotting ML ----
    # um_tree <- force.ultrametric(tree, method = "extend")
    # um_tree$root.edge <- 2
    # plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)

    # Add pie charts for ancestral states at internal nodes
    node_pies <- as.matrix(ace_results$lik.anc)

    # nodelabels(
    #   pie = node_pies,
    #   piecol = colors,
    #   cex = 0.2
    # )

    # Add pie charts for tips
    # tip_pies <- as.factor(character_data)
    # tip_pies <- to.matrix(tip_pies, levels(tip_pies))
    #
    # tiplabels(
    #   pie = tip_pies,
    #   piecol = colors,
    #   cex = 0.2
    # )

    # legend("topleft",
    #        legend = state_names,
    #        fill = colors,
    #        title = "Character States")

    # Get node number for MRCA of all Primates
    mrca <- MRCA(
      t$phy,
      t$dat |> filter(Order == "Primates") |> pull(tip.label))
    root_pie <- subset(node_pies, rownames(node_pies) %in% mrca)
    if (i == 9 | i == 10){
    root_pie <- tibble(dataset = i, rep = j, primateMRCA = mrca, DG = root_pie[,"D-G"], DP =  root_pie[, "D-P"], G = root_pie[,"G"], P = root_pie[,"P"], S = root_pie[,"S"])
    } else {
      root_pie <- tibble(dataset = i, rep = j, primateMRCA = mrca, DG = 0, DP =  0, G = root_pie[,"G"], P = root_pie[,"P"], S = root_pie[,"S"])
    }
    root_pies <- bind_rows(root_pies, root_pie)
    }
  }
}

root_pies <- rowid_to_column(root_pies)
root_pies <- root_pies |>
  mutate(tree = if_else(rowid %% 2 == 1, "Kuderna et al", "Olivier et al"))

#### Summarize ----
summary <- root_pies |>
  group_by(tree, dataset) |>
  summarise(meanDG = mean(DG), meanDP = mean(DP), meanG = mean(G), meanP = mean(P), meanS = mean(S))

summary <- pivot_longer(summary, cols = c("meanDG", "meanDP", "meanG", "meanP", "meanS"))

summary <- summary |>
  filter(dataset != 11) |>
  group_by(tree, dataset) |>
  arrange(desc(name)) |>
  mutate(prop = value) |>
  mutate(ypos = cumsum(prop)- 0.5 * prop ) |>
  rowwise() |>
  mutate(name = str_split(name, "mean")[[1]][2])

summary$display_order <- rep(c("Base", "Base + Additional"), 50)
summary$dataset <- rep(c("Olivier et al", "Olivier et al", "Kappeler & Pozzi", "Kappeler & Pozzi", "Lukas & Clutton-Brock", "Lukas & Clutton-Brock", "Müller & Thalman", "Müller & Thalman", "Müller & Thalman - Dispersed", "Müller & Thalman - Dispersed"), 10)

summary$dataset <- factor(summary$dataset, levels = c("Kappeler & Pozzi", "Olivier et al", "Lukas & Clutton-Brock", "Müller & Thalman", "Müller & Thalman - Dispersed"))

summary <- summary |>
  filter(tree == "Kuderna et al" | (tree == "Olivier et al" & display_order == "Base")) |>
  filter(value != 0) |>
  arrange(dataset, tree, display_order)

p <- ggplot(summary, aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="none") +
  scale_fill_manual(values = colors) +
  facet_wrap(vars(tree, display_order,dataset), nrow = 3, labeller = labeller(tree = c("Kuderna et al" = "", "Olivier et al" = ""))) +
  geom_label_repel(aes(y = ypos, label = name, , color = factor(name)), size = 4.5, nudge_x = 1) +
  scale_colour_manual(values = c("black", "black", "white","black", "black"))

p

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree", "summary", "p")))

Kuderna_et_al_summary <- pivot_longer(Kuderna_et_al_tree_res, cols = c("DG1", "DP1", "G1", "P1", "S1","DG2", "DP2", "G2", "P2", "S2")) |>
  filter(dataset != "Lukas & Clutton-Brock additional 2") |>
  select(dataset, phylo, name, value)

Kuderna_et_al_summary <- Kuderna_et_al_summary |>
  filter(value != 0) |>
  arrange(dataset, phylo) |>
  mutate(method = if_else(grepl("1", `name`), "ML", "SCM")) |>
  mutate(name = if_else(`name` == "DG1" | `name` == "DG2", "D-G", `name`)) |>
  mutate(name = if_else(`name` == "DP1" | `name` == "DP2", "D-P", `name`)) |>
  mutate(name = if_else(`name` == "G1" | `name` == "G2", "G", `name`)) |>
  mutate(name = if_else(`name` == "P1" | `name` == "P2", "P", `name`)) |>
  mutate(name = if_else(`name` == "S1" | `name` == "S2", "S", `name`)) |>
  group_by(dataset, method) |>
  arrange(desc(name)) |>
  mutate(prop = value) |>
  mutate(ypos = cumsum(prop)- 0.5 * prop )

Olivier_et_al_summary <- pivot_longer(Olivier_et_al_tree_res, cols = c("DG1", "DP1", "G1", "P1", "S1","DG2", "DP2", "G2", "P2", "S2")) |>
  filter(dataset != "Lukas & Clutton-Brock additional 2") |>
  select(dataset, phylo, name, value)

Olivier_et_al_summary <- Olivier_et_al_summary |>
  filter(value != 0) |>
  arrange(dataset, phylo) |>
  mutate(method = if_else(grepl("1", `name`), "ML", "SCM")) |>
  mutate(name = if_else(`name` == "DG1" | `name` == "DG2", "D-G", `name`)) |>
  mutate(name = if_else(`name` == "DP1" | `name` == "DP2", "D-P", `name`)) |>
  mutate(name = if_else(`name` == "G1" | `name` == "G2", "G", `name`)) |>
  mutate(name = if_else(`name` == "P1" | `name` == "P2", "P", `name`)) |>
  mutate(name = if_else(`name` == "S1" | `name` == "S2", "S", `name`)) |>
  group_by(dataset, method) |>
  arrange(desc(name)) |>
  mutate(prop = value) |>
  mutate(ypos = cumsum(prop)- 0.5 * prop )

colors <- c("skyblue", "red", "blue", "green", "orange")
state_names <- c("D-G", "D-P", "G", "P", "S")

p <- ggplot(Kuderna_et_al_summary, aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="none") +
  scale_fill_manual(values = colors) +
  facet_wrap(vars(method, dataset), nrow = 2) +
  geom_label_repel(aes(y = ypos, label = name, , color = factor(name)), size = 4.5, nudge_x = 1) +
  scale_colour_manual(values = c("black", "black", "white","black", "black"))

p

p <- ggplot(Olivier_et_al_summary, aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="none") +
  scale_fill_manual(values = colors) +
  facet_wrap(vars(method, dataset), nrow = 2) +
  geom_label_repel(aes(y = ypos, label = name, , color = factor(name)), size = 4.5, nudge_x = 1) +
  scale_colour_manual(values = c("black", "black", "white","black", "black"))

p
