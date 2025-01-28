library(tidyverse)
library(treedata.table)
library(ape)
library(ggtree)
library(castor)
library(phytools)
library(phangorn)
library(here)

colors <- c("skyblue", "orange", "green", "red","maroon","lavender")

# Plot Figure 2 ----
tree_file <- "Kuderna data_s4_fossil_calibrated_time_tree.nex.tree"
tree <- read.tree(tree_file)
outgroup <- c("Tupaia_belangeri", "Galeopterus_variegatus", "Mus_musculus", "Oryctolagus_cuniculus")
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(tree)
is.binary(tree)

# quick plot
plot.phylo(tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
# not tips are not aligned

d <- read_csv("Kuderna et al taxa.csv", col_names = TRUE) # adds in other taxonomic levels for taxa in Kuderna et al data set
d <- d |>
  rowwise() |>
  mutate(Genus = strsplit(label, "_")[[1]][1])

# make treedata.table object
t <- as.treedata.table(tree = tree, data = as.data.frame(d))

# d$dat holds the data...
t$dat <- t$dat |>
  rowwise() |>
  mutate(Genus = strsplit(tip.label, "_")[[1]][1])

#s <- t$dat |>
#  group_by(Genus) |>
#  summarize(count = n()) |>
#  filter(count > 1) |>
#  pull(Genus)

# t$phy holds the phylogeny
tree <- t$phy

# force the tree to be ultrametric for visualization
um_tree <- force.ultrametric(tree)
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Primate Phylogeny")
# now tips are aligned

p <- ggtree(um_tree, size = 0.3, layout = "fan") %<+% d + geom_tiplab(aes(label=label), size=2)
# plot with ggplot

# code below is to get nodes for MRCA of each Superfamily
# make dataframe for clade nodes
clades <- tibble(
  clade = c(unique(d$Superfamily)),
  node = NA
) |>
  filter(!is.na(clade) & clade != "Primates")

# find the most recent common ancestor for each clade
for (i in 1:length(clades$clade)) {
  clades$node[i] <- MRCA(
    tree,
    d |> filter(Superfamily == clades$clade[i]) |> pull(label)
  )
}

# code below is to get nodes for MRCA of each Genus
# make dataframe for clade nodes
genus_clades <- tibble(
  clade = c(unique(d$Genus)),
  node = NA
) |>
  filter(!is.na(clade))

# find the most recent common ancestor for each clade
for (i in 1:length(genus_clades$clade)) {
  genus_clades$node[i] <- MRCA(
    tree,
    d |> filter(Genus == genus_clades$clade[i]) |> pull(label)
  )
}

# plot SM Figure 2
colors <- c("red", "blue", "maroon", "skyblue", "green", "orange")

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

p

# load Character Data -----
dataG <- read_csv("Olivier et al 2024 data/dataG.csv", col_names = TRUE)
my_data <- dataG |>
  rowwise() |>
  mutate(`Genus` = str_split(`Genus_species`, "_")[[1]][1]) |>
  ungroup() |>
  mutate(main_SO = Main1) |>
  filter(sp.pop == "sp") |>
  select(Genus, Genus_species,  Common_name, superfamily, IVSO, main_SO, calculation_main_SO, Main2, foraging_style)

write_csv(my_data, "my_data.csv") # update data

outgroup <- my_data |>
  filter(superfamily == "Lorisoidea" | superfamily == "Lemuroidea") |>
  pull(Genus_species)

summary <- my_data |>
  group_by(Genus, main_SO) |>
  summarize(count = n())

# Load phylogenies from Olivier et al----
# multiple to capture uncertainty
trees = read.nexus("Olivier et al 2024 data/vert phylo.nex")

# create consensus phylogeny for robustness checks

my_tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")
is.rooted.phylo(my_tree)
my_tree <- root(my_tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(my_tree)
is.binary(my_tree)
my_tree <- multi2di(my_tree) # force tree to be binary
is.binary(my_tree)

p <- ggtree(my_tree, size = 0.1) +
  geom_tiplab(size = 2) +
  geom_text(aes(label=node), size = 2)

p

# make tree into treedata.table object
t <- as.treedata.table(tree = my_tree, data = as.data.frame(my_data))

# Filter Taxa ----
# here we can filter td for taxa of interest!

data_subset <- my_data |>
  group_by(Genus) |>
  sample_n(size = 1, replace = FALSE) |> # sample 1 species per genus
  pull(`Genus_species`)
td <- t[tip.label %in% data_subset,]

# or, use all...
td <- t

# Set up dataset ----
character_data <- td$dat$main_SO
character_data <- case_when(
  character_data == "solitary" ~ "S",
  character_data == "MF" ~ "P",
  character_data == "MFF" ~ "G",
  character_data == "FFMM" ~ "G",
  character_data == "FMM" ~ "G",
  character_data == "MF_MFF" ~ "G"
)
names(character_data) <- td$dat$tip.label
tree <- td$phy

# Convert states to a factor (required for discrete traits in phytools)
character_data <- as.factor(character_data)

# ML ancestral state reconstruction ----
ace_results <- ace(character_data, tree, type = "discrete", method = "ML")

# Define state names
state_names <- colnames(ace_results$lik.anc)

# Number the internal nodes
tree$node.label <- 1:tree$Nnode+Ntip(tree)  # Sequential labels from 1 to Nnode(tree)

# Extract the likelihoods for ancestral nodes and make dataframe
likelihoods <- as.data.frame(ace_results$lik.anc)
colnames(likelihoods) <- state_names  # Assign state names to the columns

# Convert likelihoods to data frame and add node numbers
likelihoods$node <- 1:tree$Nnode+Ntip(tree)

# Reshape the likelihoods into long format
likelihoods_long <- likelihoods |>
  pivot_longer(cols = all_of(state_names),
               names_to = "state",
               values_to = "likelihood")

# ----
# Get tree data and ensure node numbering consistency
tree_data <- fortify(tree)  # This should work now
tree_data$parent <- as.numeric(tree_data$parent)  # Ensure node numbers are numeric
likelihoods$node <- as.numeric(likelihoods$node)  # Ensure node numbers are numeric

# Merge the likelihoods with tree data by node number
tree_data <- merge(tree_data, likelihoods, by = c("parent" = "node"), all.x = TRUE)

# -----
# Force tree to be ultrametric using phytools::force_ultrametric
um_tree <- force.ultrametric(tree, method = "extend") # used by Olivier et al

p <- ggtree(um_tree,
            layout = "circular",
            size = 0.1,
            branch.length = "branch.length") +
  geom_tiplab(size = 2)

p

node_positions <- p$data |>
  filter(isTip == FALSE) |> # Only internal nodes (not tips)
  mutate(node = as.numeric(node))

# Prepare pie chart data from scratch
likelihoods_long <- likelihoods_long |>
  filter(node %in% node_positions$node)  # Only for internal nodes

likelihoods_long <- likelihoods_long |>
  left_join(node_positions, by = c("node" = "node")) |>
  mutate(x_pos = x, y_pos = y)  # Ensure the likelihoods are aligned with tree positions

pies <- nodepie(likelihoods, cols = 1:dim(likelihoods)[2] - 1)
bars <- nodebar(likelihoods, cols = 1:dim(likelihoods)[2] - 1)

# p1 <- p + geom_inset(pies, width = 0.05, height = 0.05)
# p2 <- p + geom_inset(bars, width = 0.03, height = 0.02, x = "node")

# Visualize the tree with reconstructed states using phylo.plot
# Use pie charts to display state probabilities at internal nodes
nodepies <- apply(likelihoods[, 1:dim(likelihoods)[2] - 1], 1, function(row) {
  # Each row corresponds to a node's probabilities
  row / sum(row)  # Normalize probabilities for visualization
})
colnames(nodepies) <- 1:um_tree$Nnode+Ntip(um_tree)

plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Phylogenetic Tree with Ancestral State Reconstruction")

nodelabels(pie = t(nodepies), piecol = colors, cex = 0.4)

tippies <- as.factor(character_data)
tippies <- to.matrix(tippies, levels(tippies))
tippies <- tippies[um_tree$tip.label,]

tiplabels(pie = tippies, piecol = colors, cex = 0.4)

legend("topleft", legend = state_names, fill = colors)

# SCM ancestral character state reconstruction ----

# Perform stochastic character mapping (Bayesian approach)
set.seed(123)  # For reproducibility
n_simulations <- 100  # Number of stochastic maps to generate
stochastic_maps <- make.simmap(tree, character_data, model = "ER", nsim = n_simulations, pi = "fitzjohn")

# Summarize the stochastic maps
summary_maps <- summary(stochastic_maps)

# Extract posterior probabilities for internal nodes
posterior_probs <- summary_maps$ace

# Convert observed states for tips to pie chart format
# Create a matrix for pie charts where each row represents a tip
tip_probs <- matrix(0, nrow = length(tree$tip.label), ncol = length(tree$tip.label))
state_indices <- as.numeric(character_data)  # Convert factor to numeric indices
for (i in seq_along(character_data)) {
  tip_probs[i, state_indices[i]] <- 1  # Assign 1 to the observed state
}

# Ensure posterior_probs are in a matrix format for nodes
# We need to convert the posterior probabilities into a matrix for plotting
n_nodes <- length(tree$edge[, 1])  # Number of internal nodes
node_posterior_probs <- matrix(0, nrow = n_nodes, ncol = length(state_names))

# Loop through internal nodes to assign posterior probabilities
for (i in 1:n_nodes) {
  for(j in 1: length(state_names)){
  # For each internal node, the posterior probabilities should be extracted and assigned
  node_posterior_probs[i, j] <- posterior_probs[i, j]
  }
}

# Plot the tree
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Phylogenetic Tree with Ancestral State Reconstruction")

# Add pie charts for internal nodes (posterior probabilities)
nodelabels(
  pie = node_posterior_probs,
  piecol = colors,
  cex = 0.4
)

# Add pie charts for tips (observed states)
tiplabels(
  pie = tip_probs,
  piecol = colors,
  cex = 0.4
)

# Add a legend
legend("topleft",
       legend = state_names,
       fill = colors,
       title = "Character States")

# Another ML ----

# Perform Maximum Likelihood Ancestral State Reconstruction
# Using the 'ER' model (Equal Rates model) for trait evolution
reconstructed_states <- ace(character_data, tree, type = "discrete", model = "ER")

# Extract the reconstructed states for internal nodes
reconstructed_states$lik.anc  # The reconstructed ancestral states

# Plot the tree
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Phylogenetic Tree with Ancestral State Reconstruction")

# Add pie charts for ancestral states at internal nodes
node_posterior_probs <- as.matrix(reconstructed_states$lik.anc)

nodelabels(
  pie = node_posterior_probs,
  piecol = colors,
  cex = 0.4
)

# Add pie charts for tips (observed states)
tip_probs <- matrix(0, nrow = length(tree$tip.label), ncol = length(state_names))
state_indices <- as.numeric(character_data)  # Convert factor to numeric indices
for (i in seq_along(character_data)) {
  tip_probs[i, state_indices[i]] <- 1  # Assign 1 to the observed state
}

tiplabels(
  pie = tip_probs,
  piecol = colors,
  cex = 0.4
)

# Add a legend
legend("topleft",
       legend = state_names,
       fill = colors,
       title = "Character States")



# MCMC-----

# 3. MCMC for Ancestral State Reconstruction
# Number of MCMC iterations
n_iter <- 10000
# Set initial ancestral states randomly (either 0 or 1)
initial_anc_states <- sample(c(0, 1), length(tree$node.label), replace = TRUE)

# Initialize the ancestral states vector to include both tip and internal states
# First, fill with the tip states (which are known)
anc_states <- rep(NA, length(tree$tip.label) + length(tree$node.label))
anc_states[1:length(tree$tip.label)] <- as.numeric(states)  # Tip states (0 or 1)

# 4. Set up a function for the likelihood of the model
# We will use a simple equal rates (ER) model for state transitions
log_likelihood <- function(tree, states, anc_states) {
  lik <- 0
  # Iterate over all edges in the tree
  for (i in 1:nrow(tree$edge)) {
    parent <- tree$edge[i, 1]
    child <- tree$edge[i, 2]

    # Check if the parent node is internal (>= number of tips)
    if (parent > length(tree$tip.label)) {
      # Only calculate likelihood for internal edges (not tips)
      if (anc_states[parent] != anc_states[child]) {
        lik <- lik + log(0.5)  # Transition between states (0 <-> 1) with equal probability
      }
    }
  }
  return(lik)
}

# 5. Run MCMC to sample from the posterior distribution
mcmc_samples <- matrix(NA, nrow = n_iter, ncol = length(tree$node.label))  # Store MCMC samples

# Initial state sample (ancestral)
current_states <- anc_states[(length(tree$tip.label) + 1):length(anc_states)]  # Only ancestral nodes
current_likelihood <- log_likelihood(tree, states, anc_states)

# MCMC loop
for (i in 1:n_iter) {
  # Propose a new state by randomly changing one of the ancestral states
  proposed_states <- current_states
  node_to_change <- sample(1:length(current_states), 1)
  proposed_states[node_to_change] <- sample(c(0, 1), 1)  # Randomly flip the state

  # Update the full ancestral state vector (with tips and internal nodes)
  full_proposed_states <- anc_states
  full_proposed_states[(length(tree$tip.label) + 1):length(anc_states)] <- proposed_states

  # Calculate the new likelihood
  proposed_likelihood <- log_likelihood(tree, states, full_proposed_states)

  # Accept or reject the proposed state based on the Metropolis-Hastings criterion
  acceptance_prob <- min(1, exp(proposed_likelihood - current_likelihood))

  if (runif(1) < acceptance_prob) {
    # Accept the new state
    current_states <- proposed_states
    current_likelihood <- proposed_likelihood
  }

  # Store the current state in the MCMC chain
  mcmc_samples[i, ] <- current_states
}

# 6. Analyze the MCMC samples
# The posterior distribution of states at each node
posterior_means <- apply(mcmc_samples, 2, mean)  # Mean state at each node
posterior_sd <- apply(mcmc_samples, 2, sd)  # Standard deviation (uncertainty) at each node

# 7. Plot the tree with the posterior distributions at the nodes
plot(tree, main = "Phylogenetic Tree with MCMC Ancestral State Reconstruction")

# Add pie charts for posterior state probabilities at the internal nodes
node_posterior_probs <- matrix(0, nrow = length(tree$node.label), ncol = 2)  # Two states: 0 and 1
node_posterior_probs[, 1] <- 1 - posterior_means  # State 0 probability
node_posterior_probs[, 2] <- posterior_means      # State 1 probability

# Define colors for the states
state_colors <- c("skyblue", "orange")

nodelabels(
  pie = node_posterior_probs,  # Posterior probabilities as pie chart format
  piecol = state_colors,       # Use defined colors for states
  cex = 1.2                    # Adjust size of pie charts
)

# Add pie charts for tips (observed states)
tip_probs <- matrix(0, nrow = length(tree$tip.label), ncol = 2)  # 2 states (0 and 1)
state_indices <- as.numeric(states)  # Convert factor to numeric indices
for (i in seq_along(states)) {
  tip_probs[i, state_indices[i]] <- 1  # Assign 1 to the observed state
}

tiplabels(
  pie = tip_probs,
  piecol = state_colors,
  cex = 1.2
)

# 8. Add a legend
legend("topleft", legend = c("State 0", "State 1"),
       fill = state_colors, title = "Character States")





library(ape)
is.binary(pruned_tree)
pruned_tree <- multi2di(pruned_tree)
any(pruned_tree$edge.length == 0)  # Check for zero-length branches
pruned_tree$edge.length[pruned_tree$edge.length == 0] <- 1e-8
any(is.na(d$relational_complexity))
d$relational_complexity[is.na(d$relational_complexity)] <- mean(d$relational_complexity, na.rm = TRUE)
pruned_tree <- root(pruned_tree, outgroup = "SpeciesA", resolve.root = TRUE)
reconstruction <- ace(d$relational_complexity, pruned_tree, type = "continuous", method = "ML")
reconstruction <- anc.ML(pruned_tree, d$relational_complexity)

# Step 4: View the results
print(reconstruction)

# Step 5: Plot the tree with reconstructed traits
plot(phylo_tree, main = "Tree with Ancestral Trait Reconstruction")
nodelabels(round(reconstruction$ace, 2), frame = "circle", cex = 0.8)  # Annotate nodes
tiplabels(traits, frame = "none", cex = 0.8)  # Annotate tips with observed traits

if (var(d$relational_complexity) <= 1e-8) {
  stop("Trait variance is too low for Bayesian reconstruction.")
}

set.seed(123)  # Set seed for reproducibility
reconstruction <- anc.Bayes(pruned_tree, d$relational_complexity, ngen = 10000, burnin = 2000)


# Load the necessary library
library(ape)

# Step 2: Load your species dataset
# Assuming the dataset is a vector of species names
species <- d$Species # Replace with your species list

# Step 3: Identify tips to drop
tips_to_drop <- setdiff(phylo_tree$tip.label, species)

# Step 4: Drop the tips
pruned_tree <- drop.tip(phylo_tree, tips_to_drop)

# Step 5: Visualize the pruned tree
plot(pruned_tree, main = "Pruned Phylogenetic Tree")
plot(pruned_tree, type = "fan", main = "Pruned Phylogenetic Tree")

# Step 6: Add branch length labels (optional)
nodelabels() # Adds labels to internal nodes

# Optional: Save the pruned tree to a file
write.tree(pruned_tree, file = "pruned_tree.nwk")

# Step 7: Generate the variance-covariance matrix
cov_matrix <- vcv(pruned_tree)
print(cov_matrix)


library(caper)
# Step 4: Create a comparative data object
d <- d[d$Species %in% pruned_tree$tip.label, ]
d <- d[match(pruned_tree$tip.label, d$Species), ]

comp_data <- comparative.data(pruned_tree, d, Species)

# Step 5: Fit the PGLS model
pgls_model <- pgls(relational_complexity ~ average_kinship, data = comp_data)

# Step 6: View the results
summary(pgls_model)





library(ape)
library(phytools)
# Step 1: Read the phylogenetic tree
tree_file <- "~/Desktop/newick.tree" # Replace with your tree file path
phylo_tree <- read.tree(tree_file)
phylo_tree <- phylo_tree[[1]]

## Data
library(readxl)
library(tidyverse)
d <- read_excel("~/Desktop/data.xlsx", sheet = 2, col_names = TRUE)
d[d == "NA"] <- NA
d <- d |> mutate(relational_complexity = as.numeric(relational_complexity))
d <- as.data.frame(d)

# to get primates only plus outgroup
d <- d |> filter(!is.na(`Group`) & `Species` != "Pygathrix_roxellana")

dataset <- d

matched_species <- intersect(phylo_tree$tip.label, dataset$Species)

pruned_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, matched_species))
dataset <- dataset[match(pruned_tree$tip.label, dataset$Species), ]

plot(pruned_tree, main = "Tree with Trait Values")
tiplabels(pch = 21, bg = "lightblue", cex = 1.5)  # Add tip labels
tiplabels(dataset$relational_complexity, frame = "none", adj = 1.5, cex = 0.8)  # Annotate tips with trait values

colors <- setNames(heat.colors(length(dataset$relational_complexity))[rank(dataset$relational_complexity)], dataset$Species)

# Plot tree with colored tips
plot(pruned_tree, main = "Tree with Color-Coded Traits")
tiplabels(pch = 21, bg = colors[pruned_tree$tip.label], cex = 1.5)

library(phytools)

# Check for invalid values
any(is.na(dataset$relational_complexity))  # Should return FALSE
any(is.nan(dataset$relational_complexity))  # Should return FALSE
any(is.infinite(dataset$relational_complexity))  # Should return FALSE

# Remove or handle invalid values if any
dataset <- dataset[!is.na(dataset$relational_complexity) & !is.nan(dataset$relational_complexity) & !is.infinite(dataset$relational_complexity), ]

species_to_keep <- dataset$Species
pruned_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_to_keep))
dataset <- dataset[match(pruned_tree$tip.label, dataset$Species), ]

trait_vector <- setNames(dataset$relational_complexity, dataset$Species)

is.numeric(trait_vector)  # Should return TRUE
any(is.na(trait_vector))  # Should return FALSE

cont_map <- contMap(pruned_tree, trait_vector, plot = FALSE)

# Visualize the tree with color-coded branches
plot(cont_map, ftype = "off", lwd = 3)
title("Continuous Trait Mapped onto Tree")

dataset$Species <- as.character(dataset$Species)
pruned_tree$tip.label <- tolower(pruned_tree$tip.label)
dataset$Species <- tolower(dataset$Species)
all(pruned_tree$tip.label == dataset$Species)

# Ensure species column matches tip labels exactly
dataset <- dataset[dataset$Species %in% pruned_tree$tip.label, ]

# Reorder the dataset to match tree tip order
dataset <- dataset[match(pruned_tree$tip.label, dataset$Species), ]

# Check for any leftover issues
print(dataset)

# Ensure tree is rooted
if (!is.rooted(pruned_tree)) {
  pruned_tree <- midpoint(pruned_tree)
}

# Ensure tree is fully dichotomous
if (!is.binary(pruned_tree)) {
  pruned_tree <- multi2di(pruned_tree)
}

# Fix branch lengths if any are missing or invalid
if (any(is.na(pruned_tree$edge.length)) || any(pruned_tree$edge.length <= 0)) {
  pruned_tree$edge.length[is.na(pruned_tree$edge.length) | pruned_tree$edge.length <= 0] <- 1e-8
}


# manual

aligned_species <- intersect(pruned_tree$tip.label, dataset$Species)
pruned_tree <- drop.tip(pruned_tree, setdiff(pruned_tree$tip.label, aligned_species))
dataset <- dataset[dataset$Species %in% aligned_species, ]

# Ensure correct order
dataset <- dataset[match(pruned_tree$tip.label, dataset$Species), ]

# Use dataset and pruned_tree directly in PGLS
library(nlme)

vcv_matrix <- vcv.phylo(pruned_tree)
phylo_correlation <- corBrownian(phy = pruned_tree, form = ~Species)
# Run PGLS using Brownian motion correlation structure
pgls_model <- gls(
  relational_complexity ~ average_kinship,
  correlation = phylo_correlation,
  data = dataset
)

# Summary of the PGLS model
summary(pgls_model)


# Compute the phylogenetic variance-covariance matrix
vcv_matrix <- vcv.phylo(pruned_tree)

# Calculate average kinship for each species
average_kinship <- rowMeans(vcv_matrix)
average_kinship_df <- data.frame(
  Species = names(average_kinship),
  avg_kinship = average_kinship
)

dataset$predicted <- predict(pgls_model)

# Add average kinship to dataset
dataset <- merge(dataset, average_kinship_df, by = "Species")
library(ggplot2)

# Scatter plot with regression line
ggplot(dataset |> filter(!is.na(female_infanticide)), aes(x = average_kinship, y = female_infanticide)) +
  geom_point(size = 3, color = "blue", alpha = 0.7) +  # Observed data points
  theme_minimal() +
  labs(
    x = "Average Kinship",
    y = "Trait Value",
    title = "Trait vs. Average Kinship"
  )





library(ape)
tree_file <- "genus.tree" # Replace with your tree file path

# load packages
library(ggtree)  # Bioconductor v3.10.1
library(dplyr)   # CRAN v1.1.4
library(stringr) # CRAN v1.5.1
library(geiger)  # CRAN v2.0.11

my_df <- as_tibble(phylo_tree)
my_df <- my_df %>%
  mutate(genus=str_extract(label,".*?(?=_)")) %>% # regex!
  filter(!is.na(genus)) # drops internal nodes

# process nodes for label placement
formrca <- my_df %>%
  add_count(genus, name = "numnodes") %>%
  filter(numnodes>2)

# node numbers per genus
gen_nodes <-
  formrca %>%
  group_split(genus) %>%
  purrr::map(pull, node)

# genus names
gnamesvec <- formrca %>%
  group_split(genus) %>%
  purrr::map(pull, genus) %>%
  purrr::map(unique)

# define function for getting MRCA node
getCladeNode <- function(tree, nodesvec, gname) {
  nodenum <- getMRCA(tree, tip = nodesvec)
  tibble(clade = gname, node = nodenum)
}

# table with genera and their mrca
genNodes <-
  purrr::map2_df(gen_nodes, gnamesvec,
                 ~ getCladeNode(phylo_tree, .x, .y)) # formula notation













# WORKING WITH Kuderna TREE ----

setwd('/Users/ad26693/Desktop/Active Sync/Current Writing/Annual Reviews - Evolution of Primate Social Systems')

colors <- c("skyblue", "orange", "green", "red","maroon","lavender")

tree_file <- "Kuderna data_s4_fossil_calibrated_time_tree.nex.tree" # Replace with your tree

my_tree <- read.tree(tree_file)
my_tree$tip.label
outgroup <- c(my_tree$tip.label[c(1,2,235,236)])
my_tree <- root(my_tree, outgroup = outgroup, resolve.root = TRUE)
is.rooted.phylo(my_tree)
is.binary(my_tree)

t <- as.treedata.table(tree = my_tree, data = as.data.frame(my_data))

# Filter Taxa ----

data_subset <- my_data |>
  group_by(Genus) |>
  sample_n(size = 1, replace = FALSE) |>
  pull(`Genus_species`)

td <- t[tip.label %in% data_subset,]

character_data <- td$dat$main_SO
character_data <- case_when(
  character_data == "solitary" ~ "S",
  character_data == "MF" ~ "P",
  character_data == "MFF" ~ "G",
  character_data == "FFMM" ~ "G",
  character_data == "FMM" ~ "G",
  character_data == "MF_MFF" ~ "G"
)
names(character_data) <- td$dat$tip.label
tree <- td$phy

# Convert states to a factor (required for discrete traits in phytools)
character_data <- as.factor(character_data)

#----ML ancestral state reconstruction----
ace_results <- ace(character_data, tree, type = "discrete", method = "ML")

# Define state names
state_names <- colnames(ace_results$lik.anc)

# Number the internal nodes
tree$node.label <- 1:tree$Nnode+Ntip(tree) # Sequential labels from 1 to Nnode(tree)

# Extract the likelihoods for ancestral nodes and put in dataframe
likelihoods <- as.data.frame(ace_results$lik.anc)
colnames(likelihoods) <- state_names  # Assign state names to the columns

# Add node numbers
likelihoods$node <- 1:tree$Nnode+Ntip(tree)

# Force tree to be ultrametric using phytools::force_ultrametric
um_tree <- force.ultrametric(tree, method = "extend") # used by Olivier et al

p <- ggtree(um_tree,
            layout = "circular",
            size = 0.1,
            branch.length = "branch.length"#,
            #options(ignore.negative.edge=TRUE)
) +
  geom_tiplab(size = 2)
p

# Visualize the tree with reconstructed states using phylo.plot
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "ML Ancestral State Reconstruction")

# Use pie charts to display state probabilities at internal nodes
nodepies <- apply(likelihoods[, 1:dim(likelihoods)[2] - 1], 1, function(row) {
  # Each row corresponds to a node's probabilities
  row / sum(row)  # Normalize probabilities for visualization
})
colnames(nodepies) <- 1:tree$Nnode+Ntip(tree)

# Plot node labels
# Add pie charts for internal nodes

nodelabels(pie = t(nodepies),
           piecol = colors,
           cex = 0.4)

# Make pie charts for tips
tippies <- as.factor(character_data)
tippies <- to.matrix(tippies, levels(tippies))
tippies <- tippies[um_tree$tip.label,]

# Add pie charts for tips (observed states)
tiplabels(
  pie = tippies,
  piecol = colors,
  cex = 0.4
)

# Add a legend
legend("topleft",
       legend = state_names,
       fill = colors,
       title = "Character States")

# SCM Ancestral State Reconstruction----
# Perform stochastic character mapping (Bayesian approach)
set.seed(123)  # For reproducibility
n_simulations <- 100  # Number of stochastic maps to generate
stochastic_maps <- make.simmap(tree, character_data, model = "ER", nsim = n_simulations, pi = "fitzjohn")

# Summarize the stochastic maps
summary_maps <- summary(stochastic_maps)

# Extract posterior probabilities for internal nodes
posterior_probs <- summary_maps$ace

# Convert observed states for tips to pie chart format
# Create a matrix for pie charts where each row represents a tip
tip_probs <- matrix(0, nrow = length(tree$tip.label), ncol = length(tree$tip.label))
state_indices <- as.numeric(character_data)  # Convert factor to numeric indices
for (i in seq_along(character_data)) {
  tip_probs[i, state_indices[i]] <- 1  # Assign 1 to the observed state
}

# Ensure posterior probs are in a matrix format for nodes
# We need to convert the posterior probabilities into a matrix for plotting
n_nodes <- length(tree$edge[, 1])  # Number of internal nodes
node_posterior_probs <- matrix(0, nrow = n_nodes, ncol = length(state_names))

# Loop through internal nodes to assign posterior probabilities
for (i in 1:n_nodes) {
  for(j in 1: length(state_names)){
    # For each internal node, the posterior probabilities should be extracted and assigned
    node_posterior_probs[i, j] <- posterior_probs[i, j]
  }
}

# Visualize the tree with reconstructed states using phylo.plot
plot.phylo(um_tree, type = "fan", cex = 0.6, label.offset = 4, no.margin = TRUE, main = "Stochastic Character Mapping Reconstruction")

# Add pie charts for internal nodes (posterior probabilities)
nodelabels(
  pie = node_posterior_probs,
  piecol = colors,
  cex = 0.4
)

# Add pie charts for tips (observed states)
tiplabels(
  pie = tip_probs,
  piecol = colors,
  cex = 0.4
)

# Add a legend
legend("topleft",
       legend = state_names,
       fill = colors,
       title = "Character States")

