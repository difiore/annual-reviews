# Load packages and prep workspace----
library(tidyverse)
library(treedata.table)
library(ape)
library(ggtree)
library(ggrepel)
library(phytools)
library(readxl)
library(here)

rm(list = ls()) # clear workspace
unlink("output", recursive = TRUE) # get rid of old output directory
dir.create("output") # create new output directory
# ----

# Plot SM Figure 2 - Kuderna et al phylogeny ----
tree_file <- "Kuderna data_s4_fossil_calibrated_time_tree.nex.tree"
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

p # save this figure as a PDF file 12in x 9in, open at 400dpi, and save as a PNG file 4800px x 3600px
saveRDS(p, "output/SM Figure 2.rds")
ggsave("output/SM Figure 2.pdf", width = 12, height = 9, units = "in", dpi = 400)

# clear workspace
rm(list = ls())
# ----

# Process Olivier et al data to generate common dataset ----
## load character data from Olivier et al, with additional taxonomic levels added
common_data <- read_csv("dataG for ARA.csv", col_names = TRUE)
common_data <- common_data |>
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
  select(Order, Superfamily, Species, character_data)

## save this base data file
write_csv(common_data, "output/common_data.csv")

# clear workspace
rm(list = ls())

# ----

# Data wrangling to generate datasets for analysis ----
## load in updated Olivier et al data... we keep these data for all primates other than nocturnal strepsirrhines and tarsiers
common_data <- read_csv("output/common_data.csv", col_names = TRUE)

## columns to keep
keep_columns <- c("Order", "Superfamily", "Family", "Species", "Kappeler & Pozzi 3 state","Lukas & Clutton-Brock 3 state", "Olivier et al 3 state", "Shultz et al 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 5 state", "Additional 3 state", "Additional 5 state")

## load in comparison dataset from Excel...
## we join this to base data for all primates other than nocturnal strepsirrhines and tarsiers
comparison_data <- read_excel("SM Table 2.xlsx", sheet = 1, col_names = TRUE) |>
  rowwise() |>
  mutate(`Species` = str_replace(`Species`, " ", "_")) |> # add underscore between genus and species to match to tip labels
  select(all_of(keep_columns))

combined_data <- full_join(common_data, comparison_data, by = c("Species" = "Species")) |>
  mutate(Order.y = if_else(is.na(Order.y), Order.x, Order.y)) |> # winnow to single Order column and drop extra
  mutate(Order.x = if_else(is.na(Order.x), Order.y, Order.x)) |>
  mutate(Superfamily.y = if_else(is.na(Superfamily.y), Superfamily.x, Superfamily.y)) |> # winnow to single Superfamily column and drop extra
  mutate(Superfamily.x = if_else(is.na(Superfamily.x), Superfamily.y, Superfamily.x)) |>
  select(-c("Order.y", "Superfamily.y", "Family")) |> rename(`Order` = `Order.x`, `Superfamily` = `Superfamily.x`) |> # rename single Order and Superfamily columns
  distinct() |> # keep distinct rows
  mutate(keep_index = !is.na(`character_data`) & !(Superfamily %in% c("Lemuroidea", "Lorisoidea", "Tarsioidea"))) |> # make an index to keep data from non-tarsiers and non-strepsirhines
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |> # get genus name
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
write_csv(combined_data, "output/combined_data.csv")

# clear workspace
rm(list = ls())
# ----

# Create lists of working datasets ----
## based on combined_data and additional outgroup character states
combined_data <- read_csv("output/combined_data.csv", col_names = TRUE)
d <- list() # list to hold datasets
n <- list() # list to hold dataset names

## Olivier et al 3 state
d[[1]] <- combined_data |>
  filter(!is.na(`Olivier et al 3 state`)) |>
  mutate(character_data = `Olivier et al 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
n[[1]] <- "Olivier et al 3 state"

## Olivier et al 3 state - adding in additional data for outgroup taxa
d[[2]] <- combined_data |>
  filter(!is.na(`Olivier et al 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Olivier et al 3 state`), `Olivier et al 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
n[[2]] <- "Olivier et al 3 state + outgroup"

## Kappeler & Pozzi 3 state
d[[3]] <- combined_data |>
  filter(!is.na(`Kappeler & Pozzi 3 state`)) |>
  mutate(character_data = `Kappeler & Pozzi 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
n[[3]] <- "Kappeler & Pozzi 3 state"

## Kappeler & Pozzi 3 state - adding in additional data for outgroup taxa
d[[4]] <- combined_data |>
  filter(!is.na(`Kappeler & Pozzi 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Kappeler & Pozzi 3 state`), `Kappeler & Pozzi 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
n[[4]] <- "Kappeler & Pozzi 3 state + outgroup"

## Lukas & Clutton-Brock 3 state
d[[5]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`)) |>
  mutate(character_data = `Lukas & Clutton-Brock 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
n[[5]] <- "Lukas & Clutton-Brock 3 state"

## Lukas & Clutton-Brock 3 state - adding in additional data for outgroup taxa
## note that this also replaces a few values that Lukas & Clutton-Brock had for outgroup taxa that are incorrect based on further review of the primate literature
d[[6]] <- combined_data |>
  filter(!is.na(`Lukas & Clutton-Brock 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Additional 3 state`), `Additional 3 state`, `Lukas & Clutton-Brock 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1])
n[[6]] <- "Lukas & Clutton-Brock 3 state + outgroup"

## Müller & Thalmann 3 state
d[[7]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 3 state`)) |>
  mutate(character_data = `Müller & Thalmann 3 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
n[[7]] <- "Müller & Thalmann 3 state"

## Müller & Thalmann 3 state - adding in additional data for outgroup taxa
d[[8]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 3 state`) | !is.na(`Additional 3 state`)) |>
  mutate(character_data = if_else(!is.na(`Müller & Thalmann 3 state`), `Müller & Thalmann 3 state`, `Additional 3 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "P, G", "P", `character_data`))
n[[8]] <- "Müller & Thalmann 3 state + outgroup"

## Müller & Thalmann 5 state
d[[9]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 5 state`)) |>
  mutate(character_data = `Müller & Thalmann 5 state`) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "D-P, D-G", "D-P", `character_data`)) |>
  mutate(`character_data` = if_else(`Species` == "Pongo_pygmaeus", "D-G", `character_data`))
n[[9]] <- "Müller & Thalmann 5 state"

## Müller & Thalmann 5 state - adding in additional data for outgroup taxa
d[[10]] <- combined_data |>
  filter(!is.na(`Müller & Thalmann 5 state`) | !is.na(`Additional 5 state`)) |>
  mutate(character_data = if_else(!is.na(`Müller & Thalmann 5 state`), `Müller & Thalmann 5 state`, `Additional 5 state`)) |>
  select(all_of(c("Order", "Superfamily", "Species", "character_data"))) |>
  rowwise() |>
  mutate(`Genus` = str_split(`Species`, "_")[[1]][1]) |>
  mutate(`character_data` = if_else(`character_data` == "D-P, D-G", "D-P", `character_data`)) |>
  mutate(`character_data` = if_else(`Species` == "Pongo_pygmaeus", "D-G", `character_data`))
n[[10]] <- "Müller & Thalmann 5 state + outgroup"

write_csv(d[[1]], "output/Olivier_et_al_data.csv")
write_csv(d[[2]], "output/Olivier_et_al_data_plus_outgroup.csv")

saveRDS(d, "output/datasets.rds")
saveRDS(n, "output/dataset_names.rds")

# clear workspace
rm(list = ls())
# ----

# ASR using Olivier et al phylgeny ----
## Get full tree and plot ----
## start with their original set of multiple trees to capture uncertainty
tree_file <- "vert phylo.nex"
trees = read.nexus(tree_file)
common_data <- read_csv("output/common_data.csv", col_names = TRUE)

## create consensus phylogeny using their method for generating consensus tree
tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")
is.rooted.phylo(tree)

## define outgroup for their tree, reroot, check if binary, and resolve if not
outgroup <- common_data |> filter(Superfamily %in% c("Lemuroidea", "Lorisoidea")) |> pull(Species)
is.rooted.phylo(tree)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
tree$node.label <- NULL
Olivier_et_al_tree <- tree # hold this tree!

## store tree
write.tree(tree, "output/Olivier_et_al_phylogeny.tree") # store for later

## read tree back in
tree_file <- "output/Olivier_et_al_phylogeny.tree"
tree <- read.tree(tree_file)

um_tree <- force.ultrametric(tree, method = "extend")

## test plot of full tree
p <- ggtree(um_tree,
            layout = "circular",
            size = 0.3,
            branch.length = "branch.length"
) +
  geom_tiplab(size = 2)
p

## Set up data ----
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
d <- readRDS("output/datasets.rds")
n <- readRDS("output/dataset_names.rds")

## Repeat the following for each dataset ----
### loop through all of the datasets ----
for (i in seq(1, length(d), 2)){
  # here, only loop through odd-numbered datasets...
  # i.e., no outgroups because their analysis and phylogeny did not include outgroups
  # use more colors and state_names for the Müller & Thalmann 5 state dataset
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
    pull(`tip.label`) # this is unnecessary, really, because there should be no outgroups!
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
  #### Define Data to Analyze ----
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
  # plot(fitER) # uncomment to plot transition rates
  # plot(fitARD) # uncomment to plot transition rates
  # plot(fitSYM) # uncomment to plot transition rates
  aov <- anova(fitER, fitARD, fitSYM) # compare models
  bestMk <- rownames(aov[which.min(aov$AIC),])
  Mk <- ancr(aov, type = "marginal", weighted = FALSE, tips = TRUE)
  # with weighted = false, uses the best supported model for asr
  Mk_probs <- Mk$ace
  Mk_node_pies <- Mk_probs[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), ]
  Mk_tip_pies <- Mk_probs[1:length(tree$tip.label), ]

  #### Run SCM Ancestral State Reconstruction ----
  nsim <- 100  # Number of stochastic maps to generate
  # scmER <- make.simmap(tree, character_data, model = "ER", nsim = nsim, pi = "fitzjohn")
  # scmARD <- make.simmap(tree, character_data, model = "ARD", nsim = nsim, pi = "fitzjohn")
  scmSYM <- make.simmap(tree, character_data, model = "SYM", nsim = nsim, pi = "fitzjohn")
  # SYM is best model so others commented out

  #### summarize the stochastic maps
  summary_scm <- summary(scmSYM)

  #### extract posterior probabilities and pies for internal nodes and tips
  scm_probs <- summary_scm$ace
  scm_node_pies <- scm_probs[1:tree$Nnode, ]
  scm_tip_pies <- summary_scm$tips

  #### Plotting Results ----
  ##### Mk Results ----
  # uncomment to plot Mk results
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
  # uncomment to plot SCM results, which are virtually identical to Mk results
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
    r <- tibble(dataset = n[[i]],
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
    r <- tibble(dataset = n[[i]],
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

write_csv(Olivier_et_al_tree_res, "output/Olivier_et_al_tree_res.csv")

# clear workspace
rm(list = ls())
# ----

# ASR using Kuderna et al phylogeny ----
## Get full tree and plot ----
tree_file <- "Kuderna data_s4_fossil_calibrated_time_tree.nex.tree"
tree <- read.tree(tree_file)

## replace Cephalopachus and Carlito with Tarsius to match species in Kuderna et al
tree$tip.label <- gsub("Cephalopachus", "Tarsius", tree$tip.label)
tree$tip.label <- gsub("Carlito", "Tarsius", tree$tip.label)

## define outgroup for the tree, reroot, check if binary, and resolve if not
outgroup <- c("Mus_musculus", "Oryctolagus_cuniculus", "Tupaia_belangeri", "Galeopterus_variegatus")
is.rooted.phylo(tree)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)
tree <- multi2di(tree) # force tree to be binary
is.binary(tree)
tree$node.label <- NULL
Kuderna_et_al_tree <- tree # hold this tree!

## store tree
write.tree(Kuderna_et_al_tree, "output/Kuderna_et_al_phylogeny.tree") # store for later

## read tree back in
tree_file <- "output/Kuderna_et_al_phylogeny.tree"
tree <- read.tree(tree_file)

um_tree <- force.ultrametric(tree, method = "extend")

## test plot of full tree
p <- ggtree(um_tree,
            layout = "fan",
            size = 0.3,
            branch.length = "branch.length"
) +
  geom_tiplab(size = 2)
p

## Set up data ----
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
d <- readRDS("output/datasets.rds")
n <- readRDS("output/dataset_names.rds")

## Repeat the following for each dataset ----
### loop through all of the datasets ----
for (i in seq(1, length(d))){
  # use more colors and state_names for the Müller & Thalmann 5 state datasets
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
  #### Define Data to Analyze Dataset ----
  character_data <- t$dat$`character_data`
  names(character_data) <- t$dat$tip.label
  character_data <- as.factor(character_data) # convert states to a factor (required for discrete traits in {phytools})
  tree <- t$phy

  #### Run Mk Ancestral State Reconstruction ----
  # try ER, ARD, and SYM models
  fitER <- fitMk(tree, character_data, model = "ER", pi = "fitzjohn")
  # equal rates model, fitzjohn root prior
  fitARD <- fitMk(tree, character_data, model = "ARD", pi = "fitzjohn")
  # all rates different model, fitzjohn root prior
  fitSYM <- fitMk(tree, character_data, model = "SYM", pi = "fitzjohn")
  # symmetric rates, fitzjohn root prior
  # plot(fitER) # uncomment to plot transition rates
  # plot(fitARD) # uncomment to plot transition rates
  # plot(fitSYM) # uncomment to plot transition rates
  aov <- anova(fitER, fitARD, fitSYM) # compare models
  bestMk <- rownames(aov[which.min(aov$AIC),])
  Mk <- ancr(aov, type = "marginal", weighted = FALSE, tips = TRUE)
  # with weighted = false, uses the best supported model for asr
  Mk_probs <- Mk$ace
  Mk_node_pies <- Mk_probs[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode), ]
  Mk_tip_pies <- Mk_probs[1:length(tree$tip.label), ]

  #### Run SCM Ancestral State Reconstruction ----
  nsim <- 100  # Number of stochastic maps to generate
  # scmER <- make.simmap(tree, character_data, model = "ER", nsim = nsim, pi = "fitzjohn")
  # scmARD <- make.simmap(tree, character_data, model = "ARD", nsim = nsim, pi = "fitzjohn")
  scmSYM <- make.simmap(tree, character_data, model = "SYM", nsim = nsim, pi = "fitzjohn")
  # SYM is best model so others commented out

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

  plot.phylo(um_tree,
             type = "fan",
             cex = 0.5,
             label.offset = 4,
             no.margin = TRUE,
             root.edge = TRUE)

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
  #### add legend
  legend("topright",
         inset = c(0.025),
         legend = state_names,
         cex = 0.75,
         fill = colors,
         title = "Character States")

  ##### SCM Results ----
  # uncomment to plot SCM results, which are virtually identical to Mk results
  # um_tree <- force.ultrametric(tree, method = "extend")
  # um_tree$root.edge <- 2
  #
  # plot.phylo(um_tree,
  #   type = "fan",
  #   cex = 0.5,
  #   label.offset = 4,
  #   no.margin = TRUE,
  #   root.edge = TRUE)
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

  if (i == 9 | i == 10){ # for Müller & Thalmann dataset
    r <- tibble(dataset = n[[i]],
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
  } else { # for remainder of datasets
    r <- tibble(dataset = n[[i]],
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

write_csv(Kuderna_et_al_tree_res, "output/Kuderna_et_al_tree_res.csv")

# clear workspace
rm(list = ls())
# ----

# Plot SM Figure 3 - Olivier et al data on Olivier et al tree ----
## load data
tree <- read.tree("output/Olivier_et_al_phylogeny.tree")
data <- read_csv("output/Olivier_et_al_data.csv", col_names = TRUE)

## create treedata.table to merge tree and data and drop tips that are missing in dataset and data that are missing in tree
t <- as.treedata.table(tree = tree, data = as.data.frame(data), name_column = "Species")
phy <- t$phy # get the phylogeny
is.rooted.phylo(phy) # check if it's rooted
is.binary(phy) # check if it's binary

## set outgroup to strepsirrhines
outgroup <- t$dat |>
    filter(`Superfamily` %in% c("Lemuroidea", "Lorisoidea")) |>
    pull(`tip.label`)
t$phy <- root(t$phy, outgroup = outgroup, resolve.root = TRUE) # reroot the tree
phy <- t$phy
is.rooted.phylo(phy) # check if it's rooted
is.binary(phy) # check if it's binary

## rotate some nodes for visualization
phy <- rotateNodes(phy, 219)
phy <- rotateNodes(phy, 310)
phy <- rotateNodes(phy, 325)
t$phy <- phy

### Run ASR -----
#### Set up aata ----
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
# with weighted = false, uses the best supported model for asr
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
       cex = 0.75,
       fill = colors,
       title = "Character States")

# save this figure as a PDF file 11in x 9in, open at 400dpi, and save as a PNG file 4400px x 3600px
# need to do this by hand as the plot object is a {base} R plot

# clear workspace
rm(list = ls())
# ----

# Plot SM Figure 4 - MRCA posterior probability plots ----
Olivier_et_al_tree_res <- read_csv("output/Olivier_et_al_tree_res.csv", col_names = TRUE) |>
  filter(dataset %in% c("Olivier et al 3 state", "Kappeler & Pozzi 3 state", "Lukas & Clutton-Brock 3 state", "Müller & Thalmann 3 state", "Müller & Thalmann 5 state")) |>
  filter(method1 == "Mk")

Kuderna_et_al_tree_res <- read_csv("output/Kuderna_et_al_tree_res.csv", col_names = TRUE) |>
  filter(method1 == "Mk")

combined_res <- bind_rows(Olivier_et_al_tree_res, Kuderna_et_al_tree_res)
combined_res <- pivot_longer(combined_res, cols = c("MkDG", "MkDP", "MkG", "MkP", "MkS")) |>
  select(dataset, phylo, name, value)

combined_res <- combined_res |>
  filter(value != 0) |>
  arrange(dataset) |>
  mutate(outgroup = if_else(str_detect(`dataset`, "outgroup"), "outgroup", "no outgroup")) |>
  mutate(phylo = if_else(outgroup == "outgroup", paste0(phylo, " phylogeny\n+ outgroup"), paste0(phylo, " phylogeny\n+ no outgroup"))) |>
  mutate(name = if_else(`name` == "MkDG", "D-G", `name`)) |>
  mutate(name = if_else(`name` == "MkDP", "D-P", `name`)) |>
  mutate(name = if_else(`name` == "MkG", "G", `name`)) |>
  mutate(name = if_else(`name` == "MkP", "P", `name`)) |>
  mutate(name = if_else(`name` == "MkS", "S", `name`)) |>
  group_by(phylo, outgroup, dataset) |>
  arrange(phylo, outgroup, dataset, desc(name)) |>
  mutate(prop = value) |>
  mutate(ypos = cumsum(prop)- 0.5 * prop )

combined_res$phylo <- factor(
  combined_res$phylo,
  levels = c("Olivier et al phylogeny\n+ no outgroup",
             "Kuderna et al phylogeny\n+ no outgroup",
             "Kuderna et al phylogeny\n+ outgroup"))

combined_res <- combined_res |>
  mutate(dataset = str_replace(`dataset`, " \\+ outgroup", ""))

combined_res$dataset <- factor(
  combined_res$dataset,
  levels = c("Olivier et al 3 state",
             "Kappeler & Pozzi 3 state",
             "Lukas & Clutton-Brock 3 state",
             "Müller & Thalmann 3 state",
             "Müller & Thalmann 5 state"
))

colors <- c("skyblue", "red", "blue", "green", "orange")
state_names <- c("D-G", "D-P", "G", "P", "S")

p <- ggplot(combined_res, aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_fill_manual(values = colors) +
  facet_grid(phylo ~ dataset, switch = "y") +
  geom_label_repel(aes(y = ypos, label = name, , color = factor(name)), size = 4.5, nudge_x = 1) +
  scale_colour_manual(values = c("black", "black", "grey","black", "black")) +
  labs(x = "Phylogeny", y = "Dataset") +
  scale_y_discrete(position="right", expand = c(0,0)) +
  # scale_y_continuous(sec.axis = dup_axis()) +
  theme(axis.title.x.top = element_text(size = 14, margin = margin(b=10)),
        axis.title.y.left = element_text(size = 14, margin = margin(r = 5)))

p # save this figure as a PDF file 12in x 9in, open at 400dpi, and save as a PNG file 4400px x 3600px

# save the plot object and plot
saveRDS(p, "output/SM Figure 4.rds")
ggsave("output/SM Figure 4.pdf", width = 12, height = 9, units = "in", dpi = 400)

write_csaxis.title.y.left = write_csv(combined_res, "output/combined_res.csv")

# clear workspace
rm(list=ls())
# ----
