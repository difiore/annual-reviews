# Load necessary packages
library(tidyverse)
library(treedata.table)
library(ape)
library(ggtree)
library(ggrepel)
library(phytools)
library(readxl)
library(here)

rm(list = ls()) # clear workspace

# Plot SM Figure 2 - Kuderna et al phylogeny ----
tree_file <- "Kuderna_et_al_phylogeny.tree"
tree <- read.tree(tree_file)
outgroup <- c("Tupaia_belangeri", "Galeopterus_variegatus", "Mus_musculus", "Oryctolagus_cuniculus")
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)

## quick plot
plot.phylo(tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
## note: tips are not aligned

d <- read_csv("Kuderna_et_al_taxa.csv", col_names = TRUE) # adds in other taxonomic levels for taxa in Kuderna et al data set
d <- d |>
  rowwise() |>
  mutate(Genus = strsplit(label, "_")[[1]][1])

## make treedata.table object
t <- as.treedata.table(tree = tree, data = as.data.frame(d))

## t$dat holds the data...
t$dat <- t$dat |>
  rowwise() |>
  mutate(Genus = strsplit(tip.label, "_")[[1]][1])

## t$phy holds the phylogeny
tree <- t$phy

## force the tree to be ultrametric for visualization
um_tree <- force.ultrametric(tree)
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
## now tips are aligned

## plot with ggplot
p <- ggtree(um_tree, size = 0.3, layout = "fan") %<+%
  d +
  geom_tiplab(aes(label=label), size=2)

## make a tibble to hold the node number for the MRCA of each Superfamily...
superfamilies <- tibble(
  clade = c(unique(d$Superfamily)),
  node = NA
) |>
  filter(!is.na(clade) & clade != "Primates")

## ... and then find the node number for the MRCA
for (i in 1:length(superfamilies$clade)) {
  superfamilies$node[i] <- MRCA(
    tree,
    d |> filter(Superfamily == superfamilies$clade[i]) |> pull(label)
  )
}

## make a tibble to hold the node number for the MRCA of each Genus...
genera <- tibble(
  clade = c(unique(d$Genus)),
  node = NA
) |>
  filter(!is.na(clade))

## ... and then find the node number for the MRCA
for (i in 1:length(genera$clade)) {
  genera$node[i] <- MRCA(
    tree,
    d |> filter(Genus == genera$clade[i]) |> pull(label)
  )
}

um_tree$root.edge <- 2

## final plot SM Figure 2
colors <- c("green", "blue", "orange", "red","skyblue","maroon")
p <- ggtree(um_tree, size = 0.3, layout = "fan") %<+%
  d # +
  # geom_tiplab(aes(label=label), size=2) + uncomment for tip labels
  # geom_text(aes(label=node)) # uncomment for node labels
p <- rotate(p, 417) # flip genera within Tarsioidea
p <- rotate(p, 419)
p <- rotate(p, 320) # flip clades within Catarrhini
p <- rotate(p, 321) # flip clades within Cercopithecoidea
p <- p +
  geom_highlight(data=superfamilies,
               aes(node=node,
                   fill=clade),
               alpha=0.5,
               align="left") + #,
               # show.legend = FALSE) +
  scale_fill_manual(
    values=colors,
    breaks = c("Lorisoidea", "Lemuroidea", "Tarsioidea", "Ceboidea", "Cercopithecoidea", "Hominoidea")) +
  geom_cladelab(data=genera,
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
        legend.title = element_text(size=16), # change legend title font size
        legend.text = element_text(size=12)) # change legend text font size)
p <- p + geom_rootedge()

p # save this as a high resolution PNG file 1200 x 900px

rm(list = ls())

# Process Olivier et al data to generate base dataset ----
## load character data from Olivier et al, with additional taxonomic levels added
base_data <- read_csv("dataG for ARA.csv", col_names = TRUE)
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

## save this base data file
write_csv(base_data, "base_data.csv")

rm(list = ls())

# Data wrangling to generate list of datasets for analysis ----
## load in updated Olivier et al data... we keep these base data for all primates other than nocturnal strepsirrhines and tarsiers
base_data <- read_csv("base_data.csv", col_names = TRUE)

## columns to keep
keep_columns <- c("Order", "Superfamily", "Family", "Species", "Kappeler & Pozzi 3 state","Lukas & Clutton-Brock 3 state", "Olivier et al 3 state", "Shultz et al 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 5 state", "Additional 3 state", "Additional 5 state")

## load in comparison dataset from Excel...
## we join this to base data for all primates other than nocturnal strepsirrhines and tarsiers
comparison_data <- read_excel("SM Table 2.xlsx", sheet = 1, col_names = TRUE) |>
  select(all_of(keep_columns)) |>
  mutate(`Species` = str_replace(`Species`, " ", "_")) # split scientific name
combined_data <- full_join(base_data, comparison_data, by = c("Species" = "Species")) |>
  mutate(Order.y = if_else(is.na(Order.y), Order.x, Order.y)) |> # winnow to single Order column and drop extra
  mutate(Order.x = if_else(is.na(Order.x), Order.y, Order.x)) |>
  mutate(Superfamily.y = if_else(is.na(Superfamily.y), Superfamily.x, Superfamily.y)) |> # winnow to single Superfamily column and drop extra
  mutate(Superfamily.x = if_else(is.na(Superfamily.x), Superfamily.y, Superfamily.x)) |>
  select(-c("Order.y", "Superfamily.y", "Family")) |> rename(`Order` = `Order.x`, `Superfamily` = `Superfamily.x`) |> # rename single Order and Superfamily columns
  distinct() |> # keep distinct rows
  mutate(keep_index = !is.na(`character_data`) & !(Superfamily %in% c("Lemuroidea", "Lorisoidea", "Tarsioidea"))) |> # make an index to keep data from non-tarsiers and non-strepsirhines
  rowwise() |>
  mutate(keep_index = if_else(keep_index == FALSE & Genus %in% c("Lemur", "Eulemur", "Hapalemur", "Prolemur", "Varecia", "Propithecus", "Indri"), TRUE, keep_index)) |> # add index to keep diurnal and cathermal strepsirhines
  # replace data in `character_data` with data from the appropriate column from the comparison_data table
  mutate(`Olivier et al 3 state` = if_else(keep_index == TRUE & is.na(`Olivier et al 3 state`), `character_data`, `Olivier et al 3 state`)) |>
  mutate(`Kappeler & Pozzi 3 state` = if_else(keep_index == TRUE & is.na(`Kappeler & Pozzi 3 state`), `character_data`, `Kappeler & Pozzi 3 state`)) |>
  mutate(`Lukas & Clutton-Brock 3 state` = if_else(keep_index == TRUE & is.na(`Lukas & Clutton-Brock 3 state`), `character_data`, `Lukas & Clutton-Brock 3 state`)) |>
  mutate(`Shultz et al 3 state` = if_else(keep_index == TRUE & is.na(`Shultz et al 3 state`), `character_data`, `Shultz et al 3 state`)) |>
  mutate(`Müller & Thalmann 3 state` = if_else(keep_index == TRUE & is.na(`Müller & Thalmann 3 state`), `character_data`, `Müller & Thalmann 3 state`)) |>
  mutate(`Müller & Thalmann 5 state` = if_else(keep_index == TRUE & is.na(`Müller & Thalmann 5 state`), `character_data`, `Müller & Thalmann 5 state`)) |>
  select(-c("character_data", "keep_index"))

## save this data file
write_csv(combined_data, "combined_data.csv")

rm(list = ls())

# Create lists of working datasets ----
## based on combined_data and additional outgroup character states
combined_data <- read_csv("combined_data.csv", col_names = TRUE)
d <- list() # list to hold datasets
s <- list() # list to hold study names

## Olivier et al 3 state
d[[1]] <- combined_data |>
  filter(!is.na(`Olivier et al 3 state`)) |>
  mutate(character_data = `Olivier et al 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[1]] <- "Olivier et al 3 state"

## Olivier et al 3 state - adding in additional data for outgroup taxa
d[[2]] <- combined_data |>
  filter(!is.na(`Olivier et al 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Olivier et al 3 state`), `Olivier et al 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[2]] <- "Olivier et al 3 state + outgroup"

## Kappeler & Pozzi 3 state
d[[3]] <- combined_data |>
  filter(!is.na(`Kappeler & Pozzi 3 state`)) |>
  mutate(character_data = `Kappeler & Pozzi 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[3]] <- "Kappeler & Pozzi 3 state"

## Kappeler & Pozzi 3 state - adding in additional data for outgroup taxa
d[[4]] <- combined_data |>
  filter(!is.na(`Kappeler & Pozzi 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Kappeler & Pozzi 3 state`), `Kappeler & Pozzi 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[4]] <- "Kappeler & Pozzi 3 state + outgroup"

## Lukas & Clutton-Brock 3 state
d[[5]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`)) |>
  mutate(character_data = `Lukas & Clutton-Brock 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[5]] <- "Lukas & Clutton-Brock 3 state"

## Lukas & Clutton-Brock 3 state - adding in additional data for outgroup taxa
## note that this also replaces a few values that Lukas & Clutton-Brock had for outgroup taxa that are incorrect based on further review of the primate literature
d[[6]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Additional 3 state`), `Additional 3 state`, `Lukas & Clutton-Brock 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
s[[6]] <- "Lukas & Clutton-Brock 3 state + outgroup"

## Müller & Thalmann 3 state
d[[7]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 3 state`)) |>
  mutate(character_data = `Müller & Thalmann 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[7]] <- "Müller & Thalmann 3 state"

## Müller & Thalmann 3 state - adding in additional data for outgroup taxa
d[[8]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Müller & Thalmann 3 state`), `Müller & Thalmann 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
s[[8]] <- "Müller & Thalmann 3 state + outgroup"

## Müller & Thalmann 5 state
d[[9]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 5 state`)) |>
  mutate(character_data = `Müller & Thalmann 5 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "D-P, D-G", "D-P", `character_data`)) |>
  mutate(`character_data` = if_else(`Species` == "Pongo_pygmaeus", "D-G", `character_data`))
s[[9]] <- "Müller & Thalmann 5 state"

## Müller & Thalmann 5 state - adding in additional data for outgroup taxa
d[[10]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 5 state`) | !is.na(`Additional 5 state`)) |>
  mutate(character_data = if_else(!is.na(`Müller & Thalmann 5 state`), `Müller & Thalmann 5 state`, `Additional 5 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "D-P, D-G", "D-P", `character_data`)) |>
  mutate(`character_data` = if_else(`Species` == "Pongo_pygmaeus", "D-G", `character_data`))
s[[10]] <- "Müller & Thalmann 5 state + outgroup"

rm(list = setdiff(ls(), c("d", "s")))

# ASR using Kuderna et al phylogeny ----
## Get full tree and plot ----
tree_file <- "Kuderna_et_al_phylogeny.tree"
tree <- read.tree(tree_file)

## replace Cephalopachus and Carlito with Tarsius to match species in Kuderna et al
tree$tip.label <- gsub("Cephalopachus", "Tarsius", tree$tip.label)
tree$tip.label <- gsub("Carlito", "Tarsius", tree$tip.label)

## create a vector with outgroup taxa
outgroup <- c("Mus_musculus", "Oryctolagus_cuniculus", "Tupaia_belangeri", "Galeopterus_variegatus")
is.rooted.phylo(tree)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
Kuderna_et_al_tree <- tree # hold the original tree

um_tree <- force.ultrametric(tree, method = "extend")

## test plot of full tree
p <- ggtree(um_tree,
            layout = "fan",
            size = 0.3,
            branch.length = "branch.length"
) +
  geom_tiplab(size = 2)
p

## Repeat the following for each dataset ----
### define a data structure to hold ASR results
Kuderna_et_al_tree_res <- tibble(dataset = character(),
                                 ntaxa = numeric(),
                                 ntips = numeric(),
                                 phylo = character(),
                                 method1 = character(),
                                 model1 = character(),
                                 MkDG = numeric(),
                                 MkDP = numeric(),
                                 MkG = numeric(),
                                 MkP = numeric(),
                                 MkS = numeric(),
                                 method2 = character(),
                                 model2 = character(),
                                 scmDG = numeric(),
                                 scmDP = numeric(),
                                 scmG = numeric(),
                                 scmP = numeric(),
                                 scmS = numeric())
phylo <- "Kuderna et al"

### loop through all of the datasets
for (i in seq(1, length(d))){
  # use more colors and state_names for the Müller & Thalmann 5 state datasets
  i <- 2
  if (i == 9 | i == 10) {
    colors <- c("skyblue", "red", "blue", "green")
    state_names <- c("D-G", "D-P", "G", "P") # no solitary in this dataset!
  } else {
    colors <- c("blue", "green", "orange", "red","maroon","skyblue")
    state_names <- c("G", "P", "S")
  }

  # create treedata.table to merge tree and data and drop tips that are missing in dataset and data that are missing in tree
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

  #### Run Mk Ancestral State Reconstruction ----
  # try ER, ARD, and SYM models
  fitER <- fitMk(tree, character_data, model = "ER", pi = "fitzjohn")
  # equal rates model, fitzjohn root prior
  fitARD <- fitMk(tree, character_data, model = "ARD", pi = "fitzjohn")
  # all rates different model, fitzjohn root prior
  fitSYM <- fitMk(tree, character_data, model = "SYM", pi = "fitzjohn")
  # symmetric rates, fitzjohn root prior
  # plot(fitER) # plots transition rates
  # plot(fitARD) # plots transition rates
  # plot(fitSYM) # plots transition rates
  aov <- anova(fitER, fitARD, fitSYM) # compare models
  bestMk <- rownames(aov[which.min(aov$AIC),])
  Mk <- ancr(aov, type = "marginal", weighted = FALSE, tips = TRUE)
  # with weighted = false, use the best supported model for asr
  Mk_probs <- Mk$ace
  Mk_node_pies <- Mk_probs[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), ]
  Mk_tip_pies <- Mk_probs[1:length(tree$tip.label), ]

  #### Run SCM Ancestral State Reconstruction ----
  nsim <- 100  # Number of stochastic maps to generate
  # scmER <- make.simmap(tree, character_data, model = "ER", nsim = nsim, pi = "fitzjohn")
  # scmARD <- make.simmap(tree, character_data, model = "ARD", nsim = nsim, pi = "fitzjohn")
  scmSYM <- make.simmap(tree, character_data, model = "SYM", nsim = nsim, pi = "fitzjohn")

  #### summarize the stochastic maps
  summary_scm <- summary(scmSYM)

  #### extract posterior probabilities and pies for internal nodes and tips
  scm_probs <- summary_scm$ace
  scm_node_pies <- scm_probs[1:tree$Nnode, ]
  scm_tip_pies <- summary_scm$tips

  #### Plotting Results ----
  ##### Mk Results ----
  um_tree <- force.ultrametric(tree, method = "extend")
  um_tree$root.edge <- 2

  plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)

  #### add pie charts for ancestral states at internal nodes
   nodelabels(
     pie = Mk_node_pies,
     piecol = colors,
     cex = 0.2
   )

  #### add pie charts for tips
  tiplabels(
    pie = Mk_tip_pies,
    piecol = colors,
    cex = 0.2
  )

  # legend("topleft",
  #        legend = state_names,
  #        fill = colors,
  #        title = "Character States")

  ##### SCM Results ----
  # um_tree <- force.ultrametric(tree, method = "extend")
  # um_tree$root.edge <- 2
  #
  # plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)
  #
  # #### add pie charts for ancestral states at internal nodes
  # nodelabels(
  #   pie = scm_node_pies,
  #   piecol = colors,
  #   cex = 0.2
  # )
  #
  # #### add pie charts for tips
  # tiplabels(
  #   pie = scm_tip_pies,
  #   piecol = colors,
  #   cex = 0.2
  # )
  #
  # legend("topleft",
  #        legend = state_names,
  #        fill = colors,
  #        title = "Character States")

  #### get node number for primate MRCA
  mrca <- MRCA(
    t$phy,
    t$dat |> filter(Order == "Primates") |> pull(tip.label))

  Mk_root_pie <- subset(Mk_node_pies, rownames(Mk_node_pies) %in% mrca)
  scm_root_pie <- subset(scm_node_pies, rownames(scm_node_pies) %in% mrca)

  if (i == 9 | i == 10){
    r <- tibble(dataset = s[[i]],
                phylo = phylo,
                ntaxa = nrow(t$dat),
                ntips = length(t$phy$tip.label),
                method1 = "Mk",
                model1 = bestMk,
                MkDG = Mk_root_pie[,"D-G"],
                MkDP = Mk_root_pie[,"D-P"],
                MkG = Mk_root_pie[,"G"],
                MkP = Mk_root_pie[,"P"],
                MkS = 0,
                method2 = "SCM",
                model2 = "scmSYM",
                scmDG = scm_root_pie[,"D-G"],
                scmDP = scm_root_pie[,"D-P"],
                scmG = scm_root_pie[, "G"],
                scmP = scm_root_pie[ ,"P"],
                scmS = 0)
  } else {
    r <- tibble(dataset = s[[i]],
                phylo = phylo,
                ntaxa = nrow(t$dat),
                ntips = length(t$phy$tip.label),
                method1 = "Mk",
                model1 = bestMk,
                MkDG = 0,
                MkDP = 0,
                MkG = Mk_root_pie[,"G"],
                MkP = Mk_root_pie[,"P"],
                MkS = Mk_root_pie[,"S"],
                method2 = "SCM",
                model2 = "scmSYM",
                scmDG = 0,
                scmDP = 0,
                scmG = scm_root_pie[, "G"],
                scmP = scm_root_pie[ ,"P"],
                scmS = scm_root_pie[ ,"S"])
  }

  #### store results
  Kuderna_et_al_tree_res <- bind_rows(Kuderna_et_al_tree_res, r)
}

write_csv(Kuderna_et_al_tree_res, "Kuderna_et_al_tree_res.csv")

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree")))

# ASR using Olivier et al phylgeny ----
## Get full tree and plot ----
## start with multiple trees to capture uncertainty
tree_file <- "Olivier_et_al_phylogeny.nex"
trees = read.nexus(tree_file)
base_data <- read_csv("base_data.csv", col_names = TRUE)

## create consensus phylogeny
tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")
is.rooted.phylo(tree)
outgroup <- base_data |> filter(Superfamily %in% c("Lemuroidea", "Lorisoidea")) |> pull(Species)

tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
tree$node.label <- NULL
Olivier_et_al_tree <- tree # hold the original tree

um_tree <- force.ultrametric(tree, method = "extend")

## test plot of full tree
p <- ggtree(um_tree,
            layout = "circular",
            size = 0.3,
            branch.length = "branch.length"
) +
  geom_tiplab(size = 2)
p

## Repeat the following for each dataset ----
### define a data structure to hold ASR results
Olivier_et_al_tree_res <- tibble(dataset = character(),
                                 ntaxa = numeric(),
                                 ntips = numeric(),
                                 phylo = character(),
                                 method1 = character(),
                                 model1 = character(),
                                 MkDG = numeric(),
                                 MkDP = numeric(),
                                 MkG = numeric(),
                                 MkP = numeric(),
                                 MkS = numeric(),
                                 method2 = character(),
                                 model2 = character(),
                                 scmDG = numeric(),
                                 scmDP = numeric(),
                                 scmG = numeric(),
                                 scmP = numeric(),
                                 scmS = numeric())
phylo <- "Olivier et al"

### loop through all of the datasets
for (i in seq(1, length(d), 2)){
  # use more colors and state_names for the Müller & Thalmann 5 state datasets
  if (i == 9){
    colors <- c("skyblue", "red", "blue", "green") # no solitary in this dataset!
    state_names <- c("D-G", "D-P", "G", "P")
  } else {
    colors <- c("blue", "green", "orange", "red","maroon","skyblue")
    state_names <- c("G", "P", "S")
  }

  # create treedata.table to merge tree and data and drop tips that are missing in dataset and data that are missing in tree
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

  #### Run Mk Ancestral State Reconstruction ----
  # try ER, ARD, and SYM models
  fitER <- fitMk(tree, character_data, model = "ER", pi = "fitzjohn")
  # equal rates model, fitzjohn root prior
  fitARD <- fitMk(tree, character_data, model = "ARD", pi = "fitzjohn")
  # all rates different model, fitzjohn root prior
  fitSYM <- fitMk(tree, character_data, model = "SYM", pi = "fitzjohn")
  # symmetric rates, fitzjohn root prior
  # plot(fitER) # plots transition rates
  # plot(fitARD) # plots transition rates
  # plot(fitSYM) # plots transition rates
  aov <- anova(fitER, fitARD, fitSYM) # compare models
  bestMk <- rownames(aov[which.min(aov$AIC),])
  Mk <- ancr(aov, type = "marginal", weighted = FALSE, tips = TRUE)
  # with weighted = false, use the best supported model for asr
  Mk_probs <- Mk$ace
  Mk_node_pies <- Mk_probs[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), ]
  Mk_tip_pies <- Mk_probs[1:length(tree$tip.label), ]

  #### Run SCM Ancestral State Reconstruction ----
  nsim <- 100  # Number of stochastic maps to generate
  # scmER <- make.simmap(tree, character_data, model = "ER", nsim = nsim, pi = "fitzjohn")
  # scmARD <- make.simmap(tree, character_data, model = "ARD", nsim = nsim, pi = "fitzjohn")
  scmSYM <- make.simmap(tree, character_data, model = "SYM", nsim = nsim, pi = "fitzjohn")

  #### summarize the stochastic maps
  summary_scm <- summary(scmSYM)

  #### extract posterior probabilities and pies for internal nodes and tips
  scm_probs <- summary_scm$ace
  scm_node_pies <- scm_probs[1:tree$Nnode, ]
  scm_tip_pies <- summary_scm$tips

  #### Plotting Results ----
  ##### Mk Results ----
  # um_tree <- force.ultrametric(tree, method = "extend")
  # um_tree$root.edge <- 2
  #
  # plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)
  #
  # #### add pie charts for ancestral states at internal nodes
  # nodelabels(
  #   pie = Mk_node_pies,
  #   piecol = colors,
  #   cex = 0.2
  # )
  #
  # #### add pie charts for tips
  # tiplabels(
  #   pie = Mk_tip_pies,
  #   piecol = colors,
  #   cex = 0.2
  # )
  #
  # legend("topleft",
  #        legend = state_names,
  #        fill = colors,
  #        title = "Character States")

  ##### SCM Results----
  # um_tree <- force.ultrametric(tree, method = "extend")
  # um_tree$root.edge <- 2
  #
  # plot.phylo(um_tree, type = "fan", cex = 0.5, label.offset = 4, main = "Phylogenetic Tree with Ancestral State Reconstruction", no.margin = TRUE, root.edge = TRUE)
  #
  # #### add pie charts for ancestral states at internal nodes
  # nodelabels(
  #   pie = scm_node_pies,
  #   piecol = colors,
  #   cex = 0.2
  # )
  #
  # #### add pie charts for tips
  # tiplabels(
  #   pie = scm_tip_pies,
  #   piecol = colors,
  #   cex = 0.2
  # )
  #
  # legend("topleft",
  #        legend = state_names,
  #        fill = colors,
  #        title = "Character States")

  #### get node number for primate MRCA
  mrca <- MRCA(
    t$phy,
    t$dat |> filter(Order == "Primates") |> pull(tip.label))

  Mk_root_pie <- subset(Mk_node_pies, rownames(Mk_node_pies) %in% mrca)
  scm_root_pie <- subset(scm_node_pies, rownames(scm_node_pies) %in% mrca)

  if (i == 9 | i == 10){
    r <- tibble(dataset = s[[i]],
                phylo = phylo,
                ntaxa = nrow(t$dat),
                ntips = length(t$phy$tip.label),
                method1 = "Mk",
                model1 = bestMk,
                MkDG = Mk_root_pie[,"D-G"],
                MkDP = Mk_root_pie[,"D-P"],
                MkG = Mk_root_pie[,"G"],
                MkP = Mk_root_pie[,"P"],
                MkS = 0,
                method2 = "SCM",
                model2 = "scmSYM",
                scmDG = scm_root_pie[,"D-G"],
                scmDP = scm_root_pie[,"D-P"],
                scmG = scm_root_pie[, "G"],
                scmP = scm_root_pie[ ,"P"],
                scmS = 0)
  } else {
    r <- tibble(dataset = s[[i]],
                phylo = phylo,
                ntaxa = nrow(t$dat),
                ntips = length(t$phy$tip.label),
                method1 = "Mk",
                model1 = bestMk,
                MkDG = 0,
                MkDP = 0,
                MkG = Mk_root_pie[,"G"],
                MkP = Mk_root_pie[,"P"],
                MkS = Mk_root_pie[,"S"],
                method2 = "SCM",
                model2 = "scmSYM",
                scmDG = 0,
                scmDP = 0,
                scmG = scm_root_pie[, "G"],
                scmP = scm_root_pie[ ,"P"],
                scmS = scm_root_pie[ ,"S"])
  }

  #### store results
  Olivier_et_al_tree_res <- bind_rows(Olivier_et_al_tree_res, r)
}

write_csv(Olivier_et_al_tree_res, "Olivier_et_al_tree_res")

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree")))

# Analyses of effects of taxon sampling ----
## here we can filter for taxa of interest or take random samples

root_pies <- tibble(dataset = numeric(), rep = numeric(), primateMRCA = numeric(), model = character(), DG = numeric(), DP = numeric(), G = numeric(), P = numeric(), S = numeric())

## generate 100 samples with 1 species per genus for each dataset and do ASR
for (i in 1:length(d)){
  for (j in 1:100) {
    for (k in c(Kuderna_et_al_tree, Olivier_et_al_tree)){
      # do for both trees...
      # note that we do *not* need to rename Tarsier genera again because we have held the Kuderna et al tree
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

    #### Run Mk Ancestral State Reconstruction ----
    fitSYM <- fitMk(tree, character_data, model = "SYM", pi = "fitzjohn")
    # symmetric rates, fitzjohn root prior

    Mk <- ancr(fitSYM, type = "marginal", tips = TRUE)
    Mk_probs <- Mk$ace
    Mk_node_pies <- Mk_probs[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), ]
    Mk_tip_pies <- Mk_probs[1:length(tree$tip.label), ]

    #### get node number for MRCA of all Primates
    mrca <- MRCA(
      t$phy,
      t$dat |> filter(Order == "Primates") |> pull(tip.label))
    root_pie <- subset(Mk_node_pies, rownames(Mk_node_pies) %in% mrca)
    if (i == 9 | i == 10){
    root_pie <- tibble(dataset = i, rep = j, primateMRCA = mrca, model = "fitSYM", DG = root_pie[,"D-G"], DP =  root_pie[, "D-P"], G = root_pie[,"G"], P = root_pie[,"P"], S = 0)
    } else {
      root_pie <- tibble(dataset = i, rep = j, primateMRCA = mrca, model = "fitSYM", DG = 0, DP =  0, G = root_pie[,"G"], P = root_pie[,"P"], S = root_pie[,"S"])
    }
    root_pies <- bind_rows(root_pies, root_pie)
    }
  }
}

root_pies <- rowid_to_column(root_pies)
root_pies <- root_pies |>
  mutate(tree = if_else(rowid %% 2 == 1, "Kuderna et al", "Olivier et al"))

write_csv(root_pies, "root_pies.csv")

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree", "root_pies")))

    #### Summarize ----
root_pies <- read_csv("root_pies.csv", col_names = TRUE)
summary <- root_pies |>
  group_by(tree, dataset) |>
  summarise(meanDG = mean(DG), meanDP = mean(DP), meanG = mean(G), meanP = mean(P), meanS = mean(S))

summary <- pivot_longer(summary, cols = c("meanDG", "meanDP", "meanG", "meanP", "meanS"))

summary <- summary |>
  group_by(tree, dataset) |>
  arrange(desc(name)) |>
  mutate(prop = value) |>
  mutate(ypos = cumsum(prop) - 0.5 * prop ) |>
  rowwise() |>
  mutate(name = str_split(name, "mean")[[1]][2])

summary$display_order <- rep(c("Base", "Base + Outgroup"), 50)
summary$dataset <- rep(c("Olivier et al 3 state", "Olivier et al 3 state", "Kappeler & Pozzi 3 state", "Kappeler & Pozzi 3 state", "Lukas & Clutton-Brock 3 state", "Lukas & Clutton-Brock 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 5 state", "Müller & Thalmann 5 state"), 10)

summary$dataset <- factor(summary$dataset, levels = c("Olivier et al 3 state", "Kappeler & Pozzi 3 state", "Lukas & Clutton-Brock 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 5 state"))

summary <- summary |>
  filter(tree == "Kuderna et al" | (tree == "Olivier et al" & display_order == "Base")) |>
  filter(value != 0) |>
  arrange(dataset, tree, display_order)

write_csv(summary, "summary.csv")

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree", "root_pies", "summary")))

# Lines below generate plots of ASR at Primate Root for 10 different datasets across two different phylogenies using Mk model with SYM rates where we sample randomly 1 species from each genus and average across 100 runs...

summary <- read_csv("summary.csv", col_name = TRUE)

colors <- c("skyblue", "red", "blue", "green", "orange")
state_names <- c("D-G", "D-P", "G", "P", "S")

p <- ggplot(summary, aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.position="none") +
  scale_fill_manual(values = colors) +
  facet_wrap(vars(tree, display_order,dataset), nrow = 3, labeller = labeller(tree = c("Kuderna et al" = "", "Olivier et al" = ""))) +
  geom_label_repel(aes(y = ypos, label = name, , color = factor(name)), size = 4.5, nudge_x = 1) +
  scale_colour_manual(values = c("black", "black", "white","black", "black"))

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree", "summary")))

# Plot SM Figure 3 - Olivier et al data on Olivier et al tree ----
## start with multiple trees to capture uncertainty
tree_file <- "Olivier_et_al_phylogeny.nex"
trees = read.nexus(tree_file)
base_data <- read_csv("base_data.csv", col_names = TRUE)

## create consensus phylogeny
tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")
is.rooted.phylo(tree)
outgroup <- base_data |> filter(Superfamily %in% c("Lemuroidea", "Lorisoidea")) |> pull(Species)

tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
tree$node.label <- NULL
Olivier_et_al_tree <- tree # hold the original tree

# create treedata.table to merge tree and data and drop tips that are missing in dataset and data that are missing in tree
tree <- Olivier_et_al_tree
t <- as.treedata.table(tree = tree, data = as.data.frame(d[[1]]), name_column = "Species")
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
phy <- rotateNodes(phy, 219)
phy <- rotateNodes(phy, 220)
phy <- rotateNodes(phy, 310)
phy <- rotateNodes(phy, 325)
t$phy <- phy

### Run ASR -----
#### Set Up Dataset ----
character_data <- t$dat$`character_data`
names(character_data) <- t$dat$tip.label
character_data <- as.factor(character_data) # convert states to a factor (required for discrete traits in phytools)
tree <- t$phy
tree$node.label <- 1:Nnode(tree) + length(tree$tip.label)

#### Run Mk Ancestral State Reconstruction ----
# try ER, ARD, and SYM models
fitER <- fitMk(tree, character_data, model = "ER", pi = "fitzjohn")
# equal rates model, fitzjohn root prior
fitARD <- fitMk(tree, character_data, model = "ARD", pi = "fitzjohn")
# all rates different model, fitzjohn root prior
fitSYM <- fitMk(tree, character_data, model = "SYM", pi = "fitzjohn")
aov <- anova(fitER, fitARD, fitSYM) # compare models
bestMk <- rownames(aov[which.min(aov$AIC),])
Mk <- ancr(aov, type = "marginal", weighted = FALSE, tips = TRUE)
# with weighted = false, use the best supported model for asr
Mk_probs <- Mk$ace
Mk_node_pies <- Mk_probs[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), ]
Mk_tip_pies <- Mk_probs[1:length(tree$tip.label), ]

um_tree <- force.ultrametric(tree, method = "extend")
um_tree$root.edge <- 2

colors <- c("blue", "green", "orange")
state_names <- c("G", "P", "S")
plot.phylo(um_tree,
           type = "fan",
           cex = 0.5,
           label.offset = 2,
           no.margin = TRUE,
           root.edge = TRUE,
           open.angle = 1
)

#### add pie charts for ancestral states at internal nodes
nodelabels(
  pie = Mk_node_pies,
  piecol = colors,
  cex = 0.2
)

#### add pie charts for tips
tiplabels(
  pie = Mk_tip_pies,
  piecol = colors,
  cex = 0.15
)

legend("topright",
       inset = c(0.025),
       legend = state_names,
       fill = colors,
       title = "Character States")

# save figure this as a high resolution PNG file 1200 x 900px

rm(list = ls())

# MRCA Posterior Probability Plots ----
olivier <- read_csv("Olivier_et_al_tree_res", col_names = TRUE) |>
  filter(dataset %in% c("Olivier et al 3 state", "Kappeler & Pozzi 3 state", "Lukas & Clutton-Brock 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 5 state")) |>
  filter(method1 == "Mk")

kuderna <- read_csv("Kuderna_et_al_tree_res.csv", col_names = TRUE) |>
  filter(method1 == "Mk")

res <- bind_rows(olivier, kuderna)
res <- pivot_longer(res, cols = c("MkDG", "MkDP", "MkG", "MkP", "MkS")) |>
  select(dataset, phylo, name, value)

res <- res |>
  filter(value != 0) |>
  arrange(dataset) |>
  mutate(outgroup = if_else(str_detect(`dataset`, "outgroup"), "outgroup", "no outgroup")) |>
  mutate(phylo = if_else(outgroup == "outgroup", paste0(phylo, " + outgroup"), paste0(phylo, " + no outgroup"))) |>
  mutate(name = if_else(`name` == "MkDG", "D-G", `name`)) |>
  mutate(name = if_else(`name` == "MkDP", "D-P", `name`)) |>
  mutate(name = if_else(`name` == "MkG", "G", `name`)) |>
  mutate(name = if_else(`name` == "MkP", "P", `name`)) |>
  mutate(name = if_else(`name` == "MkS", "S", `name`)) |>
  group_by(phylo, outgroup, dataset) |>
  arrange(phylo, outgroup, dataset, desc(name)) |>
  mutate(prop = value) |>
  mutate(ypos = cumsum(prop)- 0.5 * prop )

res$phylo <- factor(
  res$phylo,
  levels = c("Olivier et al + no outgroup",
             "Kuderna et al + no outgroup",
             "Kuderna et al + outgroup"))

res <- res |>
  mutate(dataset = str_replace(`dataset`, " \\+ outgroup", ""))

res$dataset <- factor(
  res$dataset,
  levels = c("Olivier et al 3 state",
             "Kappeler & Pozzi 3 state",
             "Lukas & Clutton-Brock 3 state",
             "Müller & Thalmann 3 state",
             "Müller & Thalmann 5 state"
))

colors <- c("skyblue", "red", "blue", "green", "orange")
state_names <- c("D-G", "D-P", "G", "P", "S")

p <- ggplot(res, aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_tree() +
  theme(legend.position="none") +
  scale_fill_manual(values = colors) +
  facet_grid(phylo ~ dataset, switch = "y") +
  geom_label_repel(aes(y = ypos, label = name, , color = factor(name)), size = 4.5, nudge_x = 1) +
  scale_colour_manual(values = c("black", "black", "white","black", "black"))

p

write_csv(res, "res.csv")

rm(list=setdiff(ls(), c("d", "s", "Kuderna_et_al_tree_res", "Kuderna_et_al_tree", "Olivier_et_al_tree_res", "Olivier_et_al_tree", "root_pies", "summary", "res")))

