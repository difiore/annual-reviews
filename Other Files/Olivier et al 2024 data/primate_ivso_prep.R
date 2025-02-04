#set working directory
setwd(" ")

######################################################################
#setup names for phylogeny

#load dataset to get phylogeny setup
#use R dataset (decision is arbitrary)
data = read.csv("db_R.csv")

#remove species without main SO
data = data[-which(grepl(" ",data$main_SO)), ]

#subset to population-level data
data = data[data$sp.pop=="pop",]

#names of species
spp = unique(data$Genus_species)

#load phylogenetic trees
#(100 w/ uncertainty)
library(ape)
trees = read.nexus("vert phylo.nex")

#select random tree
tree = trees[[1]]

#which species are missing in dataset
diff = setdiff(spp, tree$tip.label)
diff = data.frame(diff, common = unique(data[data$Genus_species %in% diff, "Common_name"]) )
write.csv(diff, "names_missing.csv") #species that need a manual name change

#load edited file with name conversions (manually checked)
name_change =  read.csv("vertlife tree names.csv")

#change species names to match phylogenetic tree
og_names = data$Genus_species[data$Genus_species %in% name_change$original.names]
data[data$Genus_species %in% name_change$original.names, "Genus_species"] =
  name_change$changed.names[match(og_names, name_change$original.names) ]

#check if any species lack phylogenetic data
spp = unique(data$Genus_species)
setdiff(spp, tree$tip.label) #should be 0
setdiff(tree$tip.label, spp) #should be 0


######################################################################
#final datasets for analysis

#load data
dataR = read.csv("db_R.csv")
dataG = read.csv("db_G.csv")
dataGR = read.csv("db_GR.csv")

#remove species without main SO (445 -> 223 sp)
dataR = dataR[-which(grepl(" ",dataR$main_SO)), ]
dataG = dataG[-which(grepl(" ",dataG$main_SO)), ]
dataGR = dataGR[-which(grepl(" ",dataGR$main_SO)), ]

#change species names to match phylogeny
og_namesR = dataR$Genus_species[dataR$Genus_species %in% name_change$original.names]
dataR[dataR$Genus_species %in% name_change$original.names, "Genus_species"] =
  name_change$changed.names[match(og_namesR, name_change$original.names) ]

og_namesG = dataG$Genus_species[dataG$Genus_species %in% name_change$original.names]
dataG[dataG$Genus_species %in% name_change$original.names, "Genus_species"] =
  name_change$changed.names[match(og_namesG, name_change$original.names) ]

og_namesGR = dataGR$Genus_species[dataGR$Genus_species %in% name_change$original.names]
dataGR[dataGR$Genus_species %in% name_change$original.names, "Genus_species"] =
  name_change$changed.names[match(og_namesGR, name_change$original.names) ]

######################################################################
#add new variables

#proportion of groups deviating from main SO
#higher values indicate more IVSO
dataR$IVSO_prop = 1 - dataR$calculation_main_SO
dataG$IVSO_prop = 1 - dataG$calculation_main_SO
dataGR$IVSO_prop = 1 - dataGR$calculation_main_SO

#convert to integer value to account for # of groups in
#binomial model
dataR$IVSOint = as.integer(dataR$IVSO_prop * dataR$Nbr_social_units)
dataG$IVSOint = as.integer(dataG$IVSO_prop * dataG$Nbr_social_units)
dataGR$IVSOint = as.integer(dataGR$IVSO_prop * dataGR$Nbr_social_units)

#the code below can be useful for running a categorical model
#of the most frequent SO for a given population. Species with
#two equally common SOs have an unresolved main/most frequent
#SO. The code below will have the same SO for Main 1 and Main 2
#if the main SO is resolved, otherwise it will have the two
#SOs listed. Multiple imputation can then be used to aggregate
#over this uncertainty by random sampling from the columns.
{
#organize columns for species with uncertainty in main SO
#safe to ignore warnings, they are addressed in the next step
library(tidyverse)

dataR = #separate out species with 2 main SOs
dataR %>% separate(main_SO, into = c("Main1", "Main2"), sep="_(?=[^_]+$)")
dataG = #separate out species with 2 main SOs
dataG %>% separate(main_SO, into = c("Main1", "Main2"), sep="_(?=[^_]+$)")
dataGR = #separate out species with 2 main SOs
dataGR %>% separate(main_SO, into = c("Main1", "Main2"), sep="_(?=[^_]+$)")

#replace Main2 with Main1 if only 1 main SO
dataR$Main2 = ifelse(is.na(dataR$Main2),dataR$Main1, dataR$Main2)
dataG$Main2 = ifelse(is.na(dataG$Main2),dataG$Main1, dataG$Main2)
dataGR$Main2 = ifelse(is.na(dataGR$Main2),dataGR$Main1, dataGR$Main2)
}

#add superfamily data
fam = read.csv("superfamily data.csv")
dataR$superfamily = fam$superfamily[match(dataR$Genus_species, fam$updated_name)]
dataG$superfamily = fam$superfamily[match(dataG$Genus_species, fam$updated_name)]
dataGR$superfamily = fam$superfamily[match(dataGR$Genus_species, fam$updated_name)]

#save datasets for analysis
write.csv(dataR, "dataR.csv")
write.csv(dataG, "dataG.csv")
write.csv(dataGR, "dataGR.csv")
