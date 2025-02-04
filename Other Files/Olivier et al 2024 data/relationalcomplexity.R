

### -----

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

