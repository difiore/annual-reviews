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
