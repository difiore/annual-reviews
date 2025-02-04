# SCM Ancestral State Reconstruction ----

# Perform stochastic character mapping (Bayesian approach)
set.seed(123)  # For reproducibility
n_simulations <- 100  # Number of stochastic maps to generate
stochastic_maps <- make.simmap(tree, character_data, type = "discrete", method = "ML", model = "ER", nsim = n_simulations, pi = "fitzjohn") # , Q = "mcmc")


# Plotting Olivier et al Data on Kuderna Tree ----

mammal_outgroup <- c("Tupaia_belangeri", "Galeopterus_variegatus", "Mus_musculus", "Oryctolagus_cuniculus")

primate_outgroup <-



  tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)

t <- as.treedata.table(tree = tree, data = as.data.frame(my_data))

# quick plot
plot.phylo(tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
# note tips are not aligned

# Filter Taxa ----

data_subset <- my_data |>
  group_by(Genus) |>
  sample_n(size = 1, replace = FALSE) |>
  pull(`Species`)

td <- t[tip.label %in% data_subset,]
td <- t

character_data <- td$dat$character_data
names(character_data) <- td$dat$tip.label
tree <- td$phy
