############################################################
#workspace prep
############################################################

#set directory
setwd(" ")

#load packages
library(brms)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(mice)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

#set path for cmdstan installation (optional)
cmdstanr::set_cmdstan_path("C:/ProgramData/anaconda3/envs/stan/Library/bin/cmdstan/")

#MCMC settings
n_iter <- 2000
n_warm <- 1000
n_chains <- 4

#load datasets
dataR = read.csv("dataR.csv")
dataG = read.csv("dataG.csv")
dataGR = read.csv("dataGR.csv")

#add redundant column for phylogenetic effects
dataR$phylo = dataR$Genus_species
dataG$phylo = dataG$Genus_species
dataGR$phylo = dataGR$Genus_species

#calculate % of species exhibiting MF and solitary primary SO
dataGsp = dataG[dataG$sp.pop=="sp",] #subset to species-level rows
length(unique(dataGsp[dataGsp$Main1=="MF" & dataGsp$Main2=="MF","Genus_species"])) / #% sp with primary SO being MF
  length(unique(dataGsp$Genus_species))
length(unique(dataGsp[dataGsp$Main1=="solitary" & dataGsp$Main2=="solitary","Genus_species"])) / #% sp with primary SO being solitary
  length(unique(dataGsp$Genus_species))

#number of species with uncertain IVSO
length(dataGsp[dataGsp$Nbr_social_units==1,"Genus_species"])

#number of species showing IVSO
length(unique(dataGsp[dataGsp$IVSO_prop>0,"Genus_species"])) / length(unique(dataGsp$Genus_species))

#remove redundant species-level rows (analyses are conducted on population data)
dataR = dataR[dataR$sp.pop=="pop",]
dataG = dataG[dataG$sp.pop=="pop",]
dataGR = dataGR[dataGR$sp.pop=="pop",]

#calculate % of populations exhibiting MF or solitary primary SO
nrow(dataGR[dataGR$Main1=="MF" & dataGR$Main2=="MF",]) / #% pop with primary SO being MF
  nrow(dataGR)
nrow(dataGR[dataGR$Main1=="solitary" & dataGR$Main2=="solitary",]) / #% pop with primary SO being MF
  nrow(dataGR)

#solitary foragers showing solitary as primary SO
nrow(dataGR[dataGR$foraging_style=="solitary",])
length(unique(dataGR[dataG$foraging_style=="solitary","Genus_species"]))
nrow(dataGR[dataGR$foraging_style=="solitary" & dataGR$Main1=="solitary" |
                      dataGR$foraging_style=="solitary" & dataGR$Main2=="solitary" ,]) #solitary foraging + solitary primary SO

#% of populations showing IVSO
nrow(dataGR[dataGR$IVSO_prop>0,]) / nrow(dataGR)

  #load phylogenies
#multiple to capture uncertainty
library(ape)
trees = read.nexus("vert phylo.nex")

#create consensus phylogeny for robustness checks
c_tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")

plot(c_tree) # added by AD

#data files without uncertain zeroes for IVSO robustness checks
dRr = dataR[dataR$Nbr_social_units>1,]
dGr = dataG[dataG$Nbr_social_units>1,]
dGRr = dataGR[dataGR$Nbr_social_units>1,]

#general weakly regularizing prior for all models
gen_prior =   c(eval(parse(text=
                      (paste0("c(",
                        paste0(text="prior(\"normal(0,1)\",
                        class = \"Intercept\", resp=\"SOcounts\", dpar=\"mu",
                        levels(as.factor(dataGR$Main1))[-5], "\")",collapse = ","), ")")))),
                     eval(parse(text=
                      paste0("c(",
                        paste0(text="prior(\"normal(0,1)\",
                        class = \"b\", resp=\"SOcounts\", dpar=\"mu",
                        levels(as.factor(dataGR$Main1))[-5], "\")",collapse = ","), ")"))),
                    eval(parse(text=
                      paste0("c(",
                        paste0(text="prior(\"exponential(2)\",
                        class = \"sd\", resp=\"SOcounts\", dpar=\"mu",
                        levels(as.factor(dataGR$Main1))[-5], "\")",collapse = ","), ")"))),
                    prior("normal(0,1)", class = "Intercept", resp = "IVSOint"),
                    prior("normal(0,1)", class = "b", resp = "IVSOint"),
                    prior("exponential(2)", class = "sd", resp = "IVSOint"))

############################################################
#Save image of phylogeny
############################################################

#without tips
png("phyloivso_radial.png", width = 5, height = 5, units = "in", res = 600)
plot.phylo(c_tree, type = "radial", show.tip.label = TRUE, edge.width = 2)
dev.off()

#with contour
library(phytools)
temp = aggregate(IVSO_prop ~ Genus_species, data =  dataGR, FUN = mean)
x = temp[,2]
names(x) = temp[,1]

png("phyloivso_radial_cont.png", width = 10, height = 10, units = "in", res = 600)
obj = contMap(phytools::force.ultrametric(c_tree, method="extend"),
              x, plot=FALSE, type ="fan", lwd=0.5)
n = length(obj$cols)
colfunc = colorRampPalette(c("purple4", "yellow"))
obj$cols[1:n]<-colfunc(n)
plot(obj, type ="fan", ftype = "off", lwd = 3, spread.labels = TRUE)
dev.off()

############################################################
#histograms of raw data
############################################################
#create copy of dataframe for manipulation
df1 = dataGR

#representation of superfamilies
library(ggplot2)

#organize predictors
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)

#get total counts of main (most frequent) SO
SOs = unique(df1$Main1)
t1 = table(df1$Main1)
t2 = table(df1$Main2)
valh = as.vector(NA); for(i in 1:length(t1)) { valh[i] = nrow(df1[df1$Main1==SOs[i] | df1$Main2==SOs[i],]) }
vall = as.vector(NA); for(i in 1:length(t1)) { vall[i] = valh[i] -  nrow(df1[df1$Main1==SOs[i] & df1$Main1 != df1$Main2 |
                                                                    df1$Main2==SOs[i] & df1$Main1 != df1$Main2,]) }
dfp = data.frame(cat = factor(SOs, levels = c("solitary","MF","MFF","FMM","FFMM")),
                 valh = valh, vall = vall)

colors = brewer.pal(n = 8, "Dark2")

hist1 =
ggplot(dfp, aes(x = cat, y = valh, alpha = 0.8)) +
  geom_col(color = "black", fill = "black")+
  xlab("\nPrimary social organization")+
  ylab("Number of populations\n")+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( size = "none", fill = "none", alpha = "none")+
  geom_col(inherit.aes=FALSE, data = dfp, aes(x = cat, y = vall, fill = cat, alpha = 0.8))

ggsave("figure1b_hist1.png", hist1, width = 6, height = 7)

library(ggpattern)
df1$pat = "x"
df1[which(df1$Nbr_social_units==1),"pat"] = "m"
dfp2 = data.frame(value = df1$IVSO_prop, pat = df1$pat)
hist2 =
  ggplot(dfp2, aes(x = value, pattern = pat)) +
           geom_histogram_pattern(binwidth = 0.05, fill = colors[6], color = "black",
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.99,
                   pattern_spacing = 0.01)+

  xlab("\nIntra-population variation in social organization (IVSO)")+
  ylab("Number of populations\n")+

    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( size = "none", fill = "none", alpha = "none")+
  scale_pattern_manual(values = c(m = "stripe", x = "none"))+
  guides(pattern = "none")

#combine
library(cowplot)
f1b = plot_grid(hist1, hist2, align = "v", ncol = 1)
ggsave("figure1b_hist.png", f1b, width = 5, height = 8)


############################################################
#random sample size to capture uncertainty in phylogeny
############################################################

#determine n of random samples from the phylogenetic
#tree to reduce error in estimated probabilities of SO
#and IVSO

#initial comparison is made using a naive model#####
#(without covariates) for computational efficiency

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + (1|gr(phylo, cov = A)) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~
            1 + (1|gr(phylo, cov = A)) )

#create copy of dataframe for manipulation
df1 = dataGR

#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#Values are 0 (no units observed), not missing or NA
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0

#create response variables for brms
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1)) #observation-level random effect

#estimate random effects model with 1 random phylogeny
m.re1 = brm(formula = SO_m + IVSO_m + set_rescor(FALSE),
            family = c(multinomial, binomial),
            data = df1, data2 = list(A = A),
            prior = gen_prior[gen_prior$class!="b",],
            backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
            control=list(adapt_delta=0.90, max_treedepth=12))
  #save
  saveRDS(m.re1, "m_re1.RDS")
  m.re1 = readRDS("m_re1.RDS")

#estimate random effects model with 5 random phylogenies
  #create data lists
  datal = rep(list(df1),5)
  phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
  phylol = lapply(phylol, ape::vcv, cor = TRUE)
  phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
                list(A = phylol[[3]]), list(A = phylol[[4]]),
                list(A = phylol[[5]]))

  m.re5 = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
            family = c(multinomial, binomial),
            data = datal, data2 = phylol,
            prior = gen_prior[gen_prior$class!="b",],
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
            control=list(adapt_delta=0.90, max_treedepth=12))
  #save
  saveRDS(m.re5, "m_re5.RDS")
  m.re5 = readRDS("m_re5.RDS")

#estimate random effects model with 10 random phylogenies
  #create data lists
  datal = rep(list(df1),10)
  phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),10)
  phylol = lapply(phylol, ape::vcv, cor = TRUE)
  phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
                list(A = phylol[[3]]), list(A = phylol[[4]]),
                list(A = phylol[[5]]), list(A = phylol[[6]]),
                list(A = phylol[[7]]), list(A = phylol[[8]]),
                list(A = phylol[[9]]), list(A = phylol[[10]]))

  m.re10 = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
                       family = c(multinomial, binomial),
                       data = datal, data2 = phylol,
                       prior = gen_prior[gen_prior$class!="b",],
                       warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                       control=list(adapt_delta=0.90, max_treedepth=12))
  #save
  saveRDS(m.re10, "m_re10.RDS")
  m.re10 = readRDS("m_re10.RDS")

#compare and plot phylogenetic predictions from naive model
  #get predictions
  p.re1 = fitted(m.re1, newdata = data.frame(SO_tot = 1), re_formula = NA, robust = TRUE)
  p.re5 = fitted(m.re5, newdata = data.frame(SO_tot = 1), re_formula = NA, robust = TRUE)
  p.re10 = fitted(m.re10, newdata = data.frame(SO_tot = 1), re_formula = NA, robust = TRUE)

  #estimates across N
  colnames(p.re1)[1:2] = c("median_1","mad_1") #median absolute deviation (robust SD)
  colnames(p.re5)[1:2] = c("median_5","mad_5")
  colnames(p.re10)[1:2] = c("median_10","mad_10")
  Ndf = rbind(p.re1[,1:2,1:6], p.re5[,1:2,1:6], p.re10[,1:2,1:6])
  write.csv(Ndf, "Ndf.csv")

  #difference in estimators between N1, N5, and N10 random trees
  Ndf_med = Ndf[c(1,3,5),]
  Ndf_mad = Ndf[c(2,4,6),]
   #N1 -
   round(Ndf_med[1,] - Ndf_med[2,], 3) #N5
   round(Ndf_med[1,] - Ndf_med[3,], 3) #N10
   #N5 -
   round(Ndf_med[2,] - Ndf_med[3,], 3) #N5-N10
   #Absolute delta values remain < 0.01

####
#to ensure predictions in the final model were not biased by
#using a small number of phylogenies, we also run a robustness
#check using the ecologically informed ancestral state model
#with 50 random phylogenies####

  #model formula
   SO_m = bf(SO_counts | trials(SO_tot) ~
               1 + Activity_pattern + Locom + logmean_bodysize + effort_wf +
               (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
   IVSO_m = bf(IVSOint | trials(SO_tot) ~
                 1 + Activity_pattern + Locom +  logmean_bodysize + effort_wf +
                 (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

   #data prep
   {#create copy of dataframe for manipulation
   df1 = dataGR

   #random phylo tree
   tree = trees[[sample(1:length(trees),1)]]
   A = vcv(tree, corr = TRUE)

   #sort population rows to match species-level phylogeny
   df1 = df1[order(match(df1$Genus_species, rownames(A))),]

   #prep response and predictor variables
   df1$Solitary[is.na(df1$Solitary)] = 0
   df1$MF[is.na(df1$MF)] = 0
   df1$MFF[is.na(df1$MFF)] = 0
   df1$FMM[is.na(df1$FMM)] = 0
   df1$FFMM[is.na(df1$FFMM)] = 0
   df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
   df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
   df1$obs = seq(1:nrow(df1))

   df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
   f_mean = aggregate(Nbr_papers_for_each_field_sites ~ superfamily, mean, data = df1)
   rownames(f_mean) = f_mean$superfamily
   df1$effort_bf = as.vector(scale(f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
   df1$effort_wf = as.vector(scale(df1$Nbr_papers_for_each_field_sites -
                                     f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
   df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                                 ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral",
                                        df1$Activity_pattern))
   df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
   df1$Locomotion = as.factor(df1$Locomotion)
   df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
   df1$logmean_bodysize =
     as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                    1, mean, na.rm = TRUE) )))

   datal = rep(list(df1),50)
   phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),50)
   lphy = list()
   for(i in 1:50){
   txt = paste0("phylol[[", i, "]]")
   lphy[[i]] = list(A = vcv(eval(parse(text = txt)), cor = TRUE))
   }
   }

   #loop and pull out posterior values
   mref_x502 = list()

   for(i in 1:50){
   m.ref        = brm(formula = SO_m + IVSO_m + set_rescor(FALSE),
                      family = c(multinomial, binomial), data = datal[[i]], data2 = lphy[[i]],
                      prior =  gen_prior,
                      warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                      control=list(adapt_delta=0.80, max_treedepth=10))

   #predictions
   pred = data.frame( fitted(m.ref, newdata = data.frame(effort_wf = 0, SO_tot = 1, Nbr_social_units = 1,
                                              Activity_pattern = "Nocturnal", Locom = "AR", logmean_bodysize = -2),
                             re_formula = NA, robust = FALSE, summary=FALSE))
   colnames(pred) = c("Solitary", "MF", "MFF", "FMM", "FFMM", "IVSOint")

   #probs
   pred.df = data.frame(
     median = apply(pred, 2, function(x) median(x)),
     sd = apply(pred, 2, function(x) mad(x)),
     lci = apply(pred, 2, function(x) quantile(x, c(0.05,0.95)))[1,],
     uci = apply(pred, 2, function(x) quantile(x, c(0.05,0.95)))[2,])
   pred.df = round(pred.df,2)

   mref_x502[[i]] = pred.df
   saveRDS(mref_x502, "phylopred_50.RDS")
   }

   #load results
   phylopred_50 = readRDS("phylopred_50.RDS")

   #plot SO probability as a function of sample size
   predl = do.call("rbind",phylopred_50)
   predl$SO = sub("(\\D+).*", "\\1", rownames(predl))
   predl$sim = rep(1:50,each = 6)

   #summary df
   n5 = sample(1:50,5, replace = F)
   n20 = sample(1:50,20, replace = F)
   n35 = sample(1:50,35, replace = F)

   predsum = data.frame(n = factor(rep(c("5", "20", "35", "50"), each = 6),
                                      levels=c("5", "20", "35", "50")),
               SO = rep(c("FFMM", "FMM", "Total IVSO", "MF","MFF","Solitary"), 4),
               med.est = c(
               aggregate(predl[predl$sim %in% n5,"median"], by = list(predl[predl$sim %in% n5,"SO"]),
                FUN = median)$x,
               aggregate(predl[predl$sim %in% n20,"median"], by = list(predl[predl$sim %in% n20,"SO"]),
                FUN = median)$x,
               aggregate(predl[predl$sim %in% n35,"median"], by = list(predl[predl$sim %in% n35,"SO"]),
                FUN = median)$x,
               aggregate(predl[,"median"], by = list(predl[,"SO"]), FUN = median)$x),

               sd.est = c(
               aggregate(predl[predl$sim %in% n5,"median"], by = list(predl[predl$sim %in% n5,"SO"]),
                FUN = sd)$x,
               aggregate(predl[predl$sim %in% n20,"median"], by = list(predl[predl$sim %in% n20,"SO"]),
                FUN = sd)$x,
               aggregate(predl[predl$sim %in% n35,"median"], by = list(predl[predl$sim %in% n35,"SO"]),
                FUN = sd)$x,
               aggregate(predl[,"median"], by = list(predl[,"SO"]), FUN = sd)$x),

               lci = c(
               aggregate(predl[predl$sim %in% n5,"median"], by = list(predl[predl$sim %in% n5,"SO"]),
                FUN = quantile, c(0.025,0.95))[,2][,1],
               aggregate(predl[predl$sim %in% n20,"median"], by = list(predl[predl$sim %in% n20,"SO"]),
                FUN = quantile, c(0.025,0.95))[,2][,1],
               aggregate(predl[predl$sim %in% n35,"median"], by = list(predl[predl$sim %in% n35,"SO"]),
                FUN = quantile, c(0.025,0.95))[,2][,1],
               aggregate(predl[,"median"], by = list(predl[,"SO"]), FUN = quantile, c(0.025,0.95))[,2][,1]),

               uci = c(
               aggregate(predl[predl$sim %in% n5,"median"], by = list(predl[predl$sim %in% n5,"SO"]),
                FUN = quantile, c(0.025,0.95))[,2][,2],
               aggregate(predl[predl$sim %in% n20,"median"], by = list(predl[predl$sim %in% n20,"SO"]),
                FUN = quantile, c(0.025,0.95))[,2][,2],
               aggregate(predl[predl$sim %in% n35,"median"], by = list(predl[predl$sim %in% n35,"SO"]),
                FUN = quantile, c(0.025,0.95))[,2][,2],
               aggregate(predl[,"median"], by = list(predl[,"SO"]), FUN = quantile, c(0.025,0.95))[,2][,2])
               )
  head(predsum)
  predsum[predsum$n == 5, 3:4] - predsum[predsum$n == 50, 3:4]

  #plot
  p.50pred =
  ggplot(predsum, aes(x = n, y = med.est, group = SO, color = SO)) +
  geom_point(size = 2, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = lci, ymax = uci, width = 1),
                 position = position_dodge(0.3))+
  ylab("Probability for an ancestral social unit\n")+
    xlab("Number of random phylogenies\n")+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14, angle = 35, vjust = 0.6),
          axis.title.x = element_blank(),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))

ggsave("p_50pred.png", p.50pred)


##############################################################
#Phylogenetic correlations (no ecological effects)
#############################################################

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + effort_wf + (1|superfamily) + (1|c|gr(phylo, cov = A)) + (1|c2|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 + effort_wf + (1|superfamily) + (1|c|gr(phylo, cov = A)) + (1|c2|Genus_species) + (1|obs) )

#create copy of dataframe for manipulation
df1 = dataGR

#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#prep response and predictor variables
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))
df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
f_mean = aggregate(Nbr_papers_for_each_field_sites ~ superfamily, mean, data = df1)
rownames(f_mean) = f_mean$superfamily
df1$effort_bf = as.vector(scale(f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$effort_wf = as.vector(scale(df1$Nbr_papers_for_each_field_sites - f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))

#create data lists
datal = rep(list(df1),5)
phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
phylol = lapply(phylol, ape::vcv, cor = TRUE)
phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
              list(A = phylol[[3]]), list(A = phylol[[4]]),
              list(A = phylol[[5]]))

#estimate model (running this x5 takes a long time depending on hardware)
#remove backend code line if cmdstan is not installed
m.cor = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
                        family = c(multinomial, binomial),
                        data = datal, data2 = phylol,
                        prior = c(gen_prior,prior("lkj(2)", class = "cor")),
                        backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                        warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                        control=list(adapt_delta=0.95, max_treedepth=10))
#save
saveRDS(m.cor, "m_cor.RDS")
m.cor = readRDS("m_cor.RDS")

#summarize
re = VarCorr(m.cor, summary = FALSE)
str(re$phylo$cor)

#MF
median(re$phylo$cor[,"muMF_SOcounts_Intercept","muMFF_SOcounts_Intercept"])
sum(re$phylo$cor[,"muMF_SOcounts_Intercept","muMFF_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$phylo$cor[,"muMF_SOcounts_Intercept","muFMM_SOcounts_Intercept"])
sum(re$phylo$cor[,"muMF_SOcounts_Intercept","muFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$phylo$cor[,"muMF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"])
sum(re$phylo$cor[,"muMF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$phylo$cor[,"muMF_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$phylo$cor[,"muMF_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)
#residual
median(re$Genus_species$cor[,"muMF_SOcounts_Intercept","muMFF_SOcounts_Intercept"])
sum(re$Genus_species$cor[,"muMF_SOcounts_Intercept","muMFF_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$Genus_species$cor[,"muMF_SOcounts_Intercept","muFMM_SOcounts_Intercept"])
sum(re$Genus_species$cor[,"muMF_SOcounts_Intercept","muFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$Genus_species$cor[,"muMF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"])
sum(re$Genus_species$cor[,"muMF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$Genus_species$cor[,"muMF_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$Genus_species$cor[,"muMF_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)

#MFF
median(re$phylo$cor[,"muMFF_SOcounts_Intercept","muFMM_SOcounts_Intercept"])
sum(re$phylo$cor[,"muMFF_SOcounts_Intercept","muFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$phylo$cor[,"muMFF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"])
sum(re$phylo$cor[,"muMFF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$phylo$cor[,"muMFF_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$phylo$cor[,"muMFF_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)
#residual
median(re$Genus_species$cor[,"muMFF_SOcounts_Intercept","muFMM_SOcounts_Intercept"])
sum(re$Genus_species$cor[,"muMFF_SOcounts_Intercept","muFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$Genus_species$cor[,"muMFF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"])
sum(re$Genus_species$cor[,"muMFF_SOcounts_Intercept","muFFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$Genus_species$cor[,"muMFF_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$Genus_species$cor[,"muMF_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)

#FMM
median(re$phylo$cor[,"muFMM_SOcounts_Intercept","muFFMM_SOcounts_Intercept"])
sum(re$phylo$cor[,"muFMM_SOcounts_Intercept","muFFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$phylo$cor[,"muFMM_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$phylo$cor[,"muFMM_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)
#residual
median(re$Genus_species$cor[,"muFMM_SOcounts_Intercept","muFFMM_SOcounts_Intercept"])
sum(re$Genus_species$cor[,"muFMM_SOcounts_Intercept","muFFMM_SOcounts_Intercept"] > 0)/nrow(re$phylo$cor)
median(re$Genus_species$cor[,"muFMM_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$Genus_species$cor[,"muFMM_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)

#FFMM
median(re$phylo$cor[,"muFFMM_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$phylo$cor[,"muFFMM_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)
#residual
median(re$Genus_species$cor[,"muFFMM_SOcounts_Intercept","IVSOint_Intercept"])
sum(re$Genus_species$cor[,"muFFMM_SOcounts_Intercept","IVSOint_Intercept"] > 0)/nrow(re$phylo$cor)

############################################################
#Ancestral state reconstruction (MULTINOMIAL)
############################################################

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + Activity_pattern + Locom + logmean_bodysize + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~
            1 + Activity_pattern + Locom +  logmean_bodysize + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

#for robustness check (simply change the model formula in the brms code below)
IVSO_mz = bf(IVSOint | trials(SO_tot) ~
            1 + Activity_pattern + Locom +  logmean_bodysize + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs),
            zi ~ (1|gr(phylo, cov = A)))

#create copy of dataframe for manipulation
df1 = dataGR
#df1 = dGRr #(use for robustness check without uncertain zeroes below)

#data prep
{#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#prep response and predictor variables
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))

df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
f_mean = aggregate(Nbr_papers_for_each_field_sites ~ superfamily, mean, data = df1)
rownames(f_mean) = f_mean$superfamily
df1$effort_bf = as.vector(scale(f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$effort_wf = as.vector(scale(df1$Nbr_papers_for_each_field_sites -
                                  f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                              ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral",
                                     df1$Activity_pattern))
df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
df1$Locomotion = as.factor(df1$Locomotion)
df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
df1$logmean_bodysize =
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                 1, mean, na.rm = TRUE) )))
#create data lists
datal = rep(list(df1),5)
phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
phylol = lapply(phylol, ape::vcv, cor = TRUE)
phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
              list(A = phylol[[3]]), list(A = phylol[[4]]),
              list(A = phylol[[5]]))
}

#estimate model
m.asr = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
          family = c(multinomial, binomial),
          data = datal, data2 = phylol,
          prior = gen_prior,
          backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
          warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
          control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.asr, "m_asr.RDS")
m.asr = readRDS("m_asr.RDS")

#robustness checks w/ consensus phylogeny

#zero-inflation
m.asrzi = brm(formula = SO_m + IVSO_mz + set_rescor(FALSE),
              family = c(multinomial, zero_inflated_binomial),
              data = df1, data2 = list(A = vcv(c_tree, corr=TRUE)),
              prior = gen_prior,
              backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
              warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
              control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.asrzi, "m_asrzi.RDS")
m.asrzi = readRDS("m_asrzi.RDS")

#no false zeroes
#make sure that df1 = dGRr is uncommented and run above
m.asrnz = brm(formula = SO_m + IVSO_m + set_rescor(FALSE),
                     family = c(multinomial, binomial),
                     data = df1, data2 = list(A = vcv(c_tree, corr=TRUE)),
                     prior = gen_prior,
                     backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                     warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                     control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.asrnz, "m_asrnz.RDS")
m.asrnz = readRDS("m_asrnz.RDS")

#summary of ASR predictions
pred = fitted(m.asr, newdata = data.frame(effort_wf = 0, SO_tot = 1, Nbr_social_units = 1, Activity_pattern = "Nocturnal",
                                        Locom = "AR", logmean_bodysize = -2), re_formula = NA, robust = FALSE)
N10_ASRpred = data.frame(pred = names(pred[,"Estimate",]), median = pred[,"Estimate",], sd = pred[,"Est.Error",])
N10_ASRpred[,2:3] = round(N10_ASRpred[,2:3], 2)
N10_ASRpred

pred = data.frame( fitted(m.asr, newdata = data.frame(effort_wf = 0, SO_tot = 1, Nbr_social_units = 1,
                                                      Activity_pattern = "Nocturnal", Locom = "AR", logmean_bodysize = -2),
                                                      re_formula = NA, robust = FALSE, summary=FALSE))

colnames(pred) = c("Solitary", "MF", "MFF", "FMM", "FFMM", "IVSOint")

#probs
pred.df = data.frame(
  median = apply(pred, 2, function(x) median(x)),
  sd = apply(pred, 2, function(x) mad(x)),
  lci = apply(pred, 2, function(x) quantile(x, c(0.05,0.95)))[1,],
  uci = apply(pred, 2, function(x) quantile(x, c(0.05,0.95)))[2,])
round(pred.df,3)

write.csv(round(pred.df,3),"asr_pred.csv")

#MF versus other primary SOs
diff = data.frame(pred$MF - pred)
diff.df = data.frame(
  median = apply(diff, 2, function(x) median(x)),
  sd = apply(diff, 2, function(x) mad(x)),
  lci = apply(diff, 2, function(x) quantile(x, c(0.05,0.95)))[1,],
  uci = apply(diff, 2, function(x) quantile(x, c(0.05,0.95)))[2,],
  pp = apply(diff, 2, function(x) sum(x>0)/length(x)))
round(diff.df,3)

write.csv(round(diff.df,3),"asr_diff.csv")

############################################################
#plot results
pred = data.frame( fitted(m.asr, newdata =
                          data.frame(effort_wf = 0, SO_tot = 1,
                                     Activity_pattern = "Nocturnal",
                                     Locom = "AR", logmean_bodysize = -2),
                          re_formula = NA, robust = TRUE,summary = FALSE))

colnames(pred) = c("Solitary", "MF", "MFF", "FMM", "FFMM", "Overall IVSO")

#wide to long
library(tidyr)
lnd1 = tidyr::gather(pred)
lnd1$key = factor(lnd1$key,
                  levels = c("Solitary","MF", "MFF", "FMM", "FFMM", "Overall IVSO"))

library(ggplot2)
library(tidybayes)

colors = brewer.pal(n = 8, "Dark2")
colors[[6]] = "#3ed8e6"

pred.med = data.frame(med = apply(pred, 2, median), key = factor(colnames(pred),
                  levels = c("Solitary","MF", "MFF", "FMM", "FFMM", "Overall IVSO")))

asr.primary =
  ggplot(data = lnd1, aes(x = value, color = key, fill = key, group = key))+
  geom_density(alpha = 0.25, aes(x = value, y = ..scaled..), size = 1)+
  geom_vline(data = pred.med, aes(xintercept = med, group = key), linetype = "dashed")+
  facet_wrap(. ~ key, nrow = 2, scales ="free")+
  scale_x_continuous(limits=c(0,1), expand = c(0.02,0), labels = c("0", "0.25", "0.5", "0.75", "1") )+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(values = colors[1:6]) +
  scale_fill_manual(values = colors[1:6]) +
  xlab("\nProbability for a social unit within an ancestral population")+
  ylab("Relative posterior density\n")+
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = 8),
          axis.title.y = element_text(face = "bold", size = 10),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10, face = "bold"),
          panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
          panel.background= element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(2, "lines"))+
  guides( size = "none", color = guide_legend(nrow=1), fill =
            guide_legend(nrow = 1, override.aes=list(color=NA)))


ggsave("asr_p.png", asr.primary, width = 7.5, height = 6)
ggsave("asr_p.tiff", asr.primary, width = 7.5, height = 6, compression = "lzw")

############################################################
#Compare body size effect between datasets (G and R)
############################################################

#model formula (compare body size and diet effects)
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + logmean_bodysize  +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 + logmean_bodysize +
              (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

#use consensus phylogeny (+ computational efficiency)
tree = phytools::consensus.edges(trees, method="mean.edge", if.absent="zero")
A = vcv(tree, corr = TRUE)

#create copy of dataframe for manipulation
df1 = dataG
df2 = dataR

#sort rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]
df2 = df2[order(match(df2$Genus_species, rownames(A))),]

#prep response and predictor variables
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))
df1$logmean_bodysize =
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_Galan,df1$mean_Female_Galan,
                                       df1$mean_other_Galan),1, mean, na.rm = TRUE) )))
df2$Solitary[is.na(df2$Solitary)] = 0
df2$MF[is.na(df2$MF)] = 0
df2$MFF[is.na(df2$MFF)] = 0
df2$FMM[is.na(df2$FMM)] = 0
df2$FFMM[is.na(df2$FFMM)] = 0
df2$SO_counts = with(df2, cbind(Solitary, MF, MFF, FMM, FFMM))
df2$SO_tot = with(df2, Solitary + MF + MFF + FMM + FFMM)
df2$obs = seq(1:nrow(df2))
df2$logmean_bodysize =
  as.vector(scale(log( 1 + apply(cbind(df2$mean_Male_Rowe,df2$mean_Female_Rowe,
                                       df2$mean_other_Rowe),1, mean, na.rm = TRUE) )))

#estimate model
m.G = brm(formula = SO_m + IVSO_m + set_rescor(FALSE),
          family = c(multinomial, binomial),
          data = df1, data2 = list(A = A),
          prior = gen_prior,
          backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
          warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
          control=list(adapt_delta=0.95, max_treedepth=10))

m.R = brm(formula = SO_m + IVSO_m + set_rescor(FALSE),
          family = c(multinomial, binomial),
          data = df2, data2 = list(A = A),
          prior = gen_prior,
          backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
          warmup = n_warm, iter=n_iter, chains = n_chains, inits = 0,
          control=list(adapt_delta=0.99, max_treedepth=10))

#save
saveRDS(m.G, "m_G.RDS")
saveRDS(m.R, "m_R.RDS")
m.G = readRDS("m_G.RDS")
m.R = readRDS("m_R.RDS")

#compare body size effect
post1 = as_draws_df(m.G,)
post2 = as_draws_df(m.R)

#compare distributions
diff_b = post1[ , grepl( "b_" , colnames( post1 ) ) ] -
          post2[ , grepl( "b_" , colnames( post2 ) ) ]
apply(diff_b, 2, median)
apply(diff_b, 2, quantile, c(0.05, 0.95))

############################################################
#Collective ecological effects (multivariate)
############################################################
#the model is first written in brms
#and then edited in Stan to allow for directly estimating
#variance explained in solitary (the base category)
#using a non-standard but mathematically equivalent
#parameterization with k rather than k-1 linear predictors.

#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~
            mo(Habitat_heterogenity) + Habitat_cat +
            foraging_style + Locom + Activity_pattern +
            fruitprop + folivprop + seedprop + animalprop +
            logmean_bodysize + effort_wf + (1|superfamily) +
            (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

IVSO_m = bf(IVSOint | trials(Nbr_social_units) ~ 1 +
             mo(Habitat_heterogenity) + Habitat_cat +
             foraging_style + Locom + Activity_pattern +
             fruitprop + folivprop + seedprop + animalprop +
             logmean_bodysize + effort_wf + (1|superfamily) +
             (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )

#for robustness check (simply change the model formula in the brms code below)
IVSO_mz = bf(IVSOint | trials(Nbr_social_units) ~ 1 +
             mo(Habitat_heterogenity) + Habitat_cat +
             foraging_style + Locom + Activity_pattern +
             fruitprop + folivprop + seedprop + animalprop +
             logmean_bodysize + effort_wf + (1|superfamily) +
             (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs),
            zi ~ (1|gr(phylo, cov = A)))

#create copy of dataframe for manipulation
df1 = dataGR
#df1 = dGRr (use for robustness check without uncertain zeroes)

#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#prep response and predictor variables
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))

df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
f_mean = aggregate(Nbr_papers_for_each_field_sites ~ superfamily, mean, data = df1)
rownames(f_mean) = f_mean$superfamily
df1$effort_bf = as.vector(scale(f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$effort_wf = as.vector(scale(df1$Nbr_papers_for_each_field_sites -
                                  f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                              ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral",
                                     df1$Activity_pattern))
df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
df1$Locomotion = as.factor(df1$Locomotion)
df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
df1$Locom = as.factor(df1$Locom)
df1$logmean_bodysize =
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                 1, mean, na.rm = TRUE) )))
df1$foraging_style = ifelse(df1$foraging_style == "one", "solitary",
                            ifelse(df1$foraging_style == "unknown", NA,
                                   ifelse(df1$foraging_style == " ", NA,
                                          df1$foraging_style )))
df1$foraging_style = relevel(as.factor(df1$foraging_style), ref = "solitary")
df1$fruitprop = df1$mean_Fruits_all/100
df1$folivprop = df1$mean_Leaves_all/100
df1$flowerprop = df1$mean_Flowers_all/100
df1$seedprop = df1$mean_Seeds_all/100
df1$animalprop = df1$mean_Animal_all/100
df1[df1$Habitat_cat=="","Habitat_cat"] = NA
df1$Habitat_cat = as.factor(df1$Habitat_cat)

#create data lists
datal = rep(list(df1),5)
phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
phylol = lapply(phylol, ape::vcv, cor = TRUE)
phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
              list(A = phylol[[3]]), list(A = phylol[[4]]),
              list(A = phylol[[5]]))

#generated dataset with imputation of missing values
#warnings are due to constants in df (safely ignore)
datal = lapply(datal, FUN = function(x) {
                  complete(mice(x[,c("Genus_species","phylo","superfamily","effort_wf",
                  "Solitary","MF","MFF","FMM","FFMM",
                  "Nbr_social_units", "IVSOint",
                  "Habitat_heterogenity","Habitat_cat",
                  "foraging_style","logmean_bodysize",
                  "Locom","Activity_pattern","fruitprop","folivprop",
                  "flowerprop","seedprop","animalprop")], m = 1))
  })

datal = lapply(datal, FUN = function(x) {x$SO_counts = with(x, cbind(Solitary, MF, MFF, FMM, FFMM))
                                         x$SO_tot = with(x, Solitary + MF + MFF + FMM + FFMM)
                                         x$obs = seq(1:nrow(x))
                                         return(x)})

#estimate model
m.full = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
                data = datal, data2 = phylol,
                prior = gen_prior,
                backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.full, "m_full.RDS")
m.full = readRDS("m_full.RDS")
summary(m.full, robust = TRUE)
full.fe = round(fixef(m.full, robust = TRUE, prob = c(0.05,0.95)),2)
write.csv(full.fe, "full_fe.csv")


#robustness checks w/ consensus phylogeny

#zero-inflation
m.fullzi = brm(formula = SO_m + IVSO_mz + set_rescor(FALSE), family = c(multinomial, zero_inflated_binomial),
                      data = datal[[1]], data2 = list(A = vcv(c_tree, corr=TRUE)),
                      prior = gen_prior,
                      backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                      warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                      control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.fullzi, "m_fullzi.RDS")
m.fullzi = readRDS("m_fullzi.RDS")

#no false zeroes
#make sure that df1 = dGRr is uncommented and run above
m.fullnz = brm(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
                        data = df1, data2 = list(A = vcv(c_tree, corr=TRUE)),
                        prior = gen_prior,
                        backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                        warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                        control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.fullnz, "m_fullnz.RDS")
m.fullnz = readRDS("m_fullnz.RDS")

#without imputation
datal2 = rep(list(df1),5)

m.full_nomi = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE), family = c(multinomial, binomial),
                  data = datal2, data2 = phylol,
                  prior = gen_prior,
                  backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                  warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                  control=list(adapt_delta=0.95, max_treedepth=10))

#save
saveRDS(m.full_nomi, "m_full_nomi.RDS")
m.full_nomi = readRDS("m_full_nomi.RDS")
summary(m.full_nomi)
full.fenomi = round(fixef(m.full_nomi, robust = TRUE, prob = c(0.05,0.95)),2)
write.csv(full.fenomi, "full_nomi_fe.csv")

#reparameterize in Stan

#get model code
df1.mi = datal[[1]]
txt = make_stancode(formula = SO_m + IVSO_m + set_rescor(FALSE),
              family = c(multinomial, binomial), data = df1.mi, data2 = list(A = A),
              prior = gen_prior)
write(txt,"brms_stan1.stan")

#get data files
datal.stan = datal
phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
phylol = lapply(phylol, ape::vcv, cor = TRUE)
for(i in 1:(length(phylol))){
     x = make_standata(formula = SO_m + IVSO_m + set_rescor(FALSE),
                        family = c(multinomial, binomial), data = datal.stan[[i]],
                        data2 = list(A = phylol[[i]]),
                        prior = gen_prior)

     #add extra variables for solitary variable
     x$N_21 = x$N_1
     x$N_22 = x$N_2
     x$N_23 = x$N_3
     x$N_24 = x$N_4
     x$M_21 = x$M_1
     x$M_22 = x$M_2
     x$M_23 = x$M_3
     x$M_24 = x$M_4
     x$J_21_SOcounts = x$J_1_SOcounts
     x$J_22_SOcounts = x$J_2_SOcounts
     x$J_23_SOcounts = x$J_3_SOcounts
     x$J_24_SOcounts = x$J_4_SOcounts
     x$Z_21_muSolitary_SOcounts_1 = x$Z_1_muMF_SOcounts_1
     x$Z_22_muSolitary_SOcounts_1 = x$Z_2_muMF_SOcounts_1
     x$Z_23_muSolitary_SOcounts_1 = x$Z_3_muMF_SOcounts_1
     x$Z_24_muSolitary_SOcounts_1 = x$Z_4_muMF_SOcounts_1
     x$K_muSolitary_SOcounts = x$K_muMF_SOcounts
     x$X_muSolitary_SOcounts = x$X_muMF_SOcounts
     x$Ksp_muSolitary_SOcounts = x$Ksp_muMF_SOcounts
     x$Imo_muSolitary_SOcounts = x$Imo_muMF_SOcounts
     x$Jmo_muSolitary_SOcounts = x$Jmo_muMF_SOcounts
     x$Xmo_muSolitary_SOcounts_1 = x$Xmo_muMF_SOcounts_1
     x$con_simo_muSolitary_SOcounts_1 = x$con_simo_muMF_SOcounts_1
     x$Lcov_23 = x$Lcov_3
    datal.stan[[i]] = x

    }

#run model

#load modified code
mod_r = rstan::stan_model("brms_reparam1.stan")

#estimate over 5 random trees
full.mod1 <- rstan::sampling(mod_r, data=datal.stan[[1]], init=0, iter = n_iter, warmup = n_warm,
                            chains = n_chains, cores = n_chains, control = list(adapt_delta=0.95, max_treedepth=10))
saveRDS(full.mod1, "full_mod1.RDS")

full.mod2 <- rstan::sampling(mod_r, data=datal.stan[[2]], init=0, iter = n_iter, warmup = n_warm,
                             chains = n_chains, cores = n_chains, control = list(adapt_delta=0.95, max_treedepth=10))
saveRDS(full.mod2, "full_mod2.RDS")

full.mod3 <- rstan::sampling(mod_r, data=datal.stan[[3]], init=0, iter = n_iter, warmup = n_warm,
                             chains = n_chains, cores = n_chains, control = list(adapt_delta=0.95, max_treedepth=10))
saveRDS(full.mod3, "full_mod3.RDS")

full.mod4 <- rstan::sampling(mod_r, data=datal.stan[[4]], init=0, iter = n_iter, warmup = n_warm,
                             chains = n_chains, cores = n_chains, control = list(adapt_delta=0.95, max_treedepth=10))
saveRDS(full.mod4, "full_mod4.RDS")

full.mod5 <- rstan::sampling(mod_r, data=datal.stan[[5]], init=0, iter = n_iter, warmup = n_warm,
                             chains = n_chains, cores = n_chains, control = list(adapt_delta=0.95, max_treedepth=10))
saveRDS(full.mod5, "full_mod5.RDS")

#load models
full.mod1 = readRDS("full_mod1.RDS")
full.mod2 = readRDS("full_mod2.RDS")
full.mod3 = readRDS("full_mod3.RDS")
full.mod4 = readRDS("full_mod4.RDS")
full.mod5 = readRDS("full_mod5.RDS")

#combine posteriors
comb.mod = rstan::sflist2stanfit(list(full.mod1,full.mod2,full.mod3,full.mod4,full.mod5))
saveRDS(comb.mod, "comb_modfull.RDS")

############################################################
#Plot variance explained
############################################################

#load combined posterior
comb.mod = readRDS("comb_modfull.RDS")
post = rstan::extract(comb.mod)

#variance explained by fixed effects
fe_v =  data.frame(Solitary = apply(post$muSolitary_SOcounts, 1, var) - apply(post$muSolitary_SOcounts2, 1, var),
                  MF = apply(post$muMF_SOcounts, 1, var) - apply(post$muMF_SOcounts2, 1, var),
                  MFF = apply(post$muMFF_SOcounts, 1, var) - apply(post$muMFF_SOcounts2, 1, var),
                  FMM = apply(post$muFMM_SOcounts, 1, var) - apply(post$muFMM_SOcounts2, 1, var),
                  FFMM = apply(post$muFFMM_SOcounts, 1, var) - apply(post$muFFMM_SOcounts2, 1, var),
                  IVSO = apply(post$mu_IVSOint, 1, var) - apply(post$mu_IVSOint2, 1, var))
#change negative values
fe_v[fe_v<0] = 0

species_v =  data.frame(Solitary = post$sd_21^2,
                        MF = post$sd_1^2,
                      MFF = post$sd_5^2,
                      FMM = post$sd_9^2,
                      FFMM = post$sd_13^2,
                      IVSO = post$sd_17^2)

phylo_v =  data.frame(Solitary = post$sd_23^2,
                        MF = post$sd_3^2,
                        MFF = post$sd_7^2,
                        FMM = post$sd_11^2,
                        FFMM = post$sd_15^2,
                        IVSO = post$sd_19^2)

superfamily_v =  data.frame(Solitary = post$sd_24^2,
                      MF = post$sd_4^2,
                      MFF = post$sd_8^2,
                      FMM = post$sd_12^2,
                      FFMM = post$sd_16^2,
                      IVSO = post$sd_20^2)

pop_v =  data.frame(Solitary = post$sd_22^2,
                            MF = post$sd_2^2,
                            MFF = post$sd_6^2,
                            FMM = post$sd_10^2,
                            FFMM = post$sd_14^2,
                            IVSO = post$sd_18^2)

#latent scale variance
logit_v = (pi/sqrt(3))^2

#calculate
phylo_r2 = phylo_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)

#summarize
apply(phylo_r2, 2, mean); apply(phylo_r2, 2, sd)

fe_r2 = fe_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(fe_r2, 2, mean); apply(fe_r2, 2, sd)

species_r2 = species_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(species_r2, 2, mean); apply(species_r2, 2, sd)

superf_r2 = superfamily_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(superf_r2, 2, mean); apply(superf_r2, 2, sd)

pop_r2 = logit_v / (fe_v + phylo_v + species_v + superfamily_v + logit_v)
apply(pop_r2, 2, mean); apply(pop_r2, 2, sd)

############################################################
#plot

#wide to long
library(tidyr)
ldf_phylo = tidyr::gather(phylo_r2)
ldf_fe = tidyr::gather(fe_r2)
ldf_species = tidyr::gather(species_r2)
ldf_superf = tidyr::gather(superf_r2)
ldf_pop = tidyr::gather(pop_r2)
ldf_phylo$type = rep("Phylogeny",nrow(ldf_phylo))
ldf_fe$type = rep("Ecological predictors",nrow(ldf_fe))
ldf_species$type = rep("Species residual",nrow(ldf_species))
ldf_superf$type = rep("Superfamily residual",nrow(ldf_superf))
ldf_pop$type = rep("Population residual",nrow(ldf_pop))

ldf = rbind(ldf_phylo, ldf_fe, ldf_species, ldf_superf, ldf_pop)
ldf$key[ldf$key=="IVSO"] = "Overall IVSO"
ldf$key = factor(ldf$key, levels = c("Solitary", "MF","MFF" , "FMM", "FFMM", "Overall IVSO"))
ldf$type = factor(ldf$type, levels = c("Phylogeny", "Ecological predictors",  "Population residual",
                                       "Species residual", "Superfamily residual"))
library(ggplot2)
library(tidybayes)

cc <- scales::seq_gradient_pal("lightblue", "darkblue")(seq(0,1,length.out=5))

R2.plot =
  ggplot(ldf, aes(x = value, y = type, color = type, fill = type))+
  coord_cartesian(xlim=c(0,1))+
  stat_pointinterval(.width = c(0.90), size = 10, position = position_dodge(width = 0.5)) +
  facet_wrap(. ~ key, nrow=1)+
  scale_color_manual(values = cc)+
  scale_fill_manual(values = cc)+
  scale_y_discrete(limits=rev(levels(ldf$type)))+
  scale_x_continuous(expand=c(0,0),breaks=c(0,0.5,1), labels =c("0","0.5","1"))+
  xlab(" \nProportion of variance explained")+
  theme(plot.title =element_text(size=12, hjust=0.5),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=10, face = "bold"),
        axis.title.y=element_blank(),
        axis.text.x=element_text(size=10),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA),
        axis.line = element_line(size = 1),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 12),
        panel.spacing = unit(1.3, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1,
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="top",
        plot.margin = unit(c(0.1,0.4,0.7,0.05), "cm"))

#save
save_plot("variance_explained.png", R2.plot, base_height = 3.5, base_width = 12)
saveRDS(R2.plot, "r2_plot.RDS")
R2.plot = readRDS("r2_plot.RDS")

############################################################
#Specific ecological effects (univariate)
############################################################
#each section fits a model for a single predictor
#and then creates a plot for the total effect

#create copy of dataframe for manipulation
df1 = dataGR

#random phylo tree
tree = trees[[sample(1:length(trees),1)]]
A = vcv(tree, corr = TRUE)

#sort population rows to match species-level phylogeny
df1 = df1[order(match(df1$Genus_species, rownames(A))),]

#prep response and predictor variables
df1$Solitary[is.na(df1$Solitary)] = 0
df1$MF[is.na(df1$MF)] = 0
df1$MFF[is.na(df1$MFF)] = 0
df1$FMM[is.na(df1$FMM)] = 0
df1$FFMM[is.na(df1$FFMM)] = 0
df1$SO_counts = with(df1, cbind(Solitary, MF, MFF, FMM, FFMM))
df1$SO_tot = with(df1, Solitary + MF + MFF + FMM + FFMM)
df1$obs = seq(1:nrow(df1))

df1$effort = as.vector(scale(log(df1$Nbr_papers_for_each_field_sites)))
f_mean = aggregate(Nbr_papers_for_each_field_sites ~ superfamily, mean, data = df1)
rownames(f_mean) = f_mean$superfamily
df1$effort_bf = as.vector(scale(f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$effort_wf = as.vector(scale(df1$Nbr_papers_for_each_field_sites -
                                  f_mean[df1$superfamily,"Nbr_papers_for_each_field_sites"] ))
df1$Activity_pattern = ifelse(df1$Activity_pattern == "Cathemeral_Diurnal_Nocturnal", "Cathemeral",
                              ifelse(df1$Activity_pattern == "Cathemeral_Diurnal", "Cathemeral",
                                     df1$Activity_pattern))
df1$Activity_pattern = relevel(as.factor(df1$Activity_pattern), ref = "Nocturnal")
df1$Locomotion = as.factor(df1$Locomotion)
df1$Locom = (ifelse(df1$Locomotion == " ", NA, as.character(df1$Locomotion)))
df1$Locom = as.factor(df1$Locom)
df1$logmean_bodysize =
  as.vector(scale(log( 1 + apply(cbind(df1$mean_Male_all,df1$mean_Female_all,df1$mean_other_all),
                                 1, mean, na.rm = TRUE) )))
df1$foraging_style = ifelse(df1$foraging_style == "one", "solitary",
                            ifelse(df1$foraging_style == "unknown", NA,
                                   ifelse(df1$foraging_style == " ", NA,
                                          df1$foraging_style )))
df1$foraging_style = relevel(as.factor(df1$foraging_style), ref = "solitary")
df1$fruitprop = df1$mean_Fruits_all/100
df1$folivprop = df1$mean_Leaves_all/100
df1$flowerprop = df1$mean_Flowers_all/100
df1$seedprop = df1$mean_Seeds_all/100
df1$animalprop = df1$mean_Animal_all/100
df1[df1$Habitat_cat=="","Habitat_cat"] = NA
df1$Habitat_cat = as.factor(df1$Habitat_cat)

#create data lists
datal = rep(list(df1),5)
phylol = rep(list(A = trees[[sample(1:length(trees),1)]]),5)
phylol = lapply(phylol, ape::vcv, cor = TRUE)
phylol = list(list(A = phylol[[1]]), list(A = phylol[[2]]),
              list(A = phylol[[3]]), list(A = phylol[[4]]),
              list(A = phylol[[5]]))


############################################################
#research effort
#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~
            1 + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))

m.eff = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
            family = c(multinomial, binomial), data = datal, data2 = phylol,
            prior = gen_prior,
            backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
            control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.eff, "m_eff.RDS")
m.eff = readRDS("m_eff.RDS")
summary(m.eff, robust = TRUE)
fixef(m.eff, robust = TRUE)


############################################################
#activity
#model formula
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + Activity_pattern + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs) )
IVSO_m = bf(IVSOint | trials(SO_tot) ~
            1 + Activity_pattern + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))

m.act = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
            family = c(multinomial, binomial), data = datal, data2 = phylol,
            prior = gen_prior,
            backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
            control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.act, "m_act.RDS")
m.act = readRDS("m_act.RDS")
summary(m.act, robust = TRUE)
fixef(m.act, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.act, resp = "SOcounts", effect = "Activity_pattern",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_Activity_pattern:cats__`$data

c2  = plot(conditional_effects(m.act, resp = "IVSOint", effect = "Activity_pattern",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_Activity_pattern$data

df2$cats__ = "Overall IVSO"
df2$effect1__ = df2$Activity_pattern
df2$effect2__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")
colors[[6]] = "#3ed8e6"

df3$response = factor(df3$effect1__, levels = c("Solitary", "MF", "MFF", "FMM", "FFMM","Overall IVSO"))

act.p1 =
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df3$lower__, ymax = df3$upper__, width = 0,
                position=position_dodge(width=0.5))+
  xlab("\nActivity pattern")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

#calculate comparisons of interest
newdf = data.frame(Activity_pattern = unique(df1$Activity_pattern), SO_tot = 1, effort_wf = 0)
pred = fitted(m.act, newdata = newdf, re_formula = NA, summary = FALSE, scale = "response")
comp1 = apply(pred, 3, FUN = function(x) x[,1] - x[,2] ) #nocturnal - diurnal
comp2 = apply(pred, 3, FUN = function(x) x[,3] - x[,2] ) #cathemeral - diurnal
comp3 = apply(pred, 3, FUN = function(x) x[,1] - x[,3] ) #nocturnal - diurnal

act_comp =
  data.frame(
    outcome = names(pred[1,1,]),
    diff_noc_diu = paste0(round(apply(comp1,2,median),2), "(",
                          apply(apply(comp1, 2, quantile, c(0.05, 0.95)),2,
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_cat_diu = paste0(round(apply(comp2,2,median),2), "(",
                          apply(apply(comp2, 2, quantile, c(0.05, 0.95)),2,
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_noc_cat = paste0(round(apply(comp3,2,median),2), "(",
                          apply(apply(comp3, 2, quantile, c(0.05, 0.95)),2,
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")") )

write.csv(act_comp,"act_comp.csv")
act_comp

############################################################
#locomotion
SO_m = bf(SO_counts | trials(SO_tot) ~
              1 + Locom + effort_wf +
             (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 + Locom + effort_wf +
              (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species)  + (1|obs))

m.loc = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
            family = c(multinomial, binomial), data = datal, data2 = phylol,
            prior = gen_prior,
            backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
            warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
            control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.loc, "m_loc.RDS")
m.loc = readRDS("m_loc.RDS")
summary(m.loc, robust = TRUE)
fixef(m.loc, robust = TRUE)

############################################################
c1  = plot(conditional_effects(m.loc, resp = "SOcounts", effect = "Locom",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_Locom:cats__`$data

c2  = plot(conditional_effects(m.loc, resp = "IVSOint", effect = "Locom",
                               conditions = data.frame(SO_tot=1),
                               robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_Locom$data

df2$cats__ = "Overall IVSO"
df2$effect1__ = df2$Locom
df2$effect2__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")
colors[[6]] = "#3ed8e6"

df3$response = factor(df3$effect1__, levels = c("Solitary", "MF", "MFF", "FMM", "FFMM","Overall IVSO"))


loc.p1 =
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_x_discrete(labels = c("Arboreal", "Both", "Terrestrial"))+
  scale_y_continuous(lim=c(0,1))+
  geom_point(position=position_dodge(width=0.5), size = 2)+
  geom_errorbar(ymin = df3$lower__, ymax = df3$upper__, width = 0,
                position=position_dodge(width=0.5))+
  xlab("Locomotion")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

#calculate comparisons of interest
newdf = data.frame(Locom = unique(df1$Locom), SO_tot = 1, effort_wf = 0)
pred = fitted(m.loc, newdata = newdf, re_formula = NA, summary = FALSE, scale = "response")
comp1 = apply(pred, 3, FUN = function(x) x[,1] - x[,2] ) #AR - BOTH
comp2 = apply(pred, 3, FUN = function(x) x[,3] - x[,4] ) #BOTH - T
comp3 = apply(pred, 3, FUN = function(x) x[,1] - x[,4] ) #AR - T

loc_comp =
  data.frame(
    outcome = names(pred[1,1,]),
    diff_AR_BOTH = paste0(round(apply(comp1,2,median),2), "(",
                          apply(apply(comp1, 2, quantile, c(0.05, 0.95)),2,
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                          ")"),
    diff_BOTH_T = paste0(round(apply(comp2,2,median),2), "(",
                         apply(apply(comp2, 2, quantile, c(0.05, 0.95)),2,
                               FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                         ")"),
    diff_AR_T = paste0(round(apply(comp3,2,median),2), "(",
                       apply(apply(comp3, 2, quantile, c(0.05, 0.95)),2,
                             FUN = function(x) paste(round(x[1],2), ",", round(x[2],2))),
                       ")") )

write.csv(loc_comp,"loc_comp.csv")
loc_comp

############################################################
#body size
SO_m = bf(SO_counts | trials(SO_tot) ~
            1 + logmean_bodysize + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 + logmean_bodysize + effort_wf +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.bs = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
           family = c(multinomial, binomial), data = datal, data2 = phylol,
           prior = gen_prior,
           backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
           warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
           control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.bs, "m_bs.RDS")
m.bs = readRDS("m_bs.RDS")
summary(m.bs, robust = TRUE)
fixef(m.bs, robust = TRUE)

############################################################
#MainSO
c1  = plot(conditional_effects(m.bs, resp = "SOcounts", effect = "logmean_bodysize",
                               conditions = data.frame(SO_tot = 1),
                               categorical = TRUE, robust = TRUE, prob = 0.9), plot = FALSE)
df = c1$`SOcounts.SOcounts_logmean_bodysize:cats__`$data

c2  = plot(conditional_effects(m.bs, resp = "IVSOint", effect = "logmean_bodysize",
                               conditions = data.frame(SO_tot=1), robust = TRUE, prob = 0.9), plot = FALSE)
df2 = c2$IVSOint.IVSOint_logmean_bodysize$data

#df$response = "x"
df2$effect2__ = "Overall IVSO"
df2$cats__ = "Overall IVSO"

df3 = rbind(df,df2)

colors = brewer.pal(n = 8, "Dark2")
colors[[6]] = "#3ed8e6"

bs.p1 =
  ggplot(df3, aes(x = effect1__, y = estimate__, group = cats__, color = cats__))+
  scale_y_continuous(lim=c(0,1))+
  facet_wrap(.~cats__, nrow = 1)+
  geom_line(size = 1)+
  geom_ribbon(aes(ymin = df3$lower__, ymax = df3$upper__, fill = cats__), color = NA, alpha = 0.2)+
  xlab("\n Log mean body size")+
  ylab("Probability\n")+
  scale_color_manual(values = colors[1:6])+
  scale_fill_manual(values = colors[1:6])+
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(face = "bold", size = 10),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10, face = "bold"),
        panel.border=element_rect(fill=NA,color="black", size=1, linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))+
  guides(colour = guide_legend(nrow = 1))

#calculate comparisons of interest
newdf = data.frame(logmean_bodysize = c(-1,0,1), SO_tot = 1, effort_wf=0)
pred = fitted(m.bs, newdata = newdf, re_formula = NA, summary = FALSE, scale = "response")
comp1 = apply(pred, 3, FUN = function(x) x[,1] - x[,2] ) #low - avg
comp2 = apply(pred, 3, FUN = function(x) x[,1] - x[,3] ) #low - high
comp3 = apply(pred, 3, FUN = function(x) x[,2] - x[,3] ) #avg - high

bs_comp =
  data.frame(
    outcome = names(pred[1,1,]),
    diff_low_avg = paste0(round(apply(comp1,2,median),2), "(",
                          apply(apply(comp1, 2, quantile, c(0.05, 0.95)),2,
                                FUN = function(x) paste(round(x[1],2), ",", round(x[2],3))),
                          ")"),
    diff_low_high = paste0(round(apply(comp2,2,median),2), "(",
                         apply(apply(comp2, 2, quantile, c(0.05, 0.95)),2,
                               FUN = function(x) paste(round(x[1],2), ",", round(x[2],3))),
                         ")"),
    diff_avg_high = paste0(round(apply(comp3,2,median),2), "(",
                           apply(apply(comp3, 2, quantile, c(0.05, 0.95)),2,
                                 FUN = function(x) paste(round(x[1],2), ",", round(x[2],3))),
                           ")"))

write.csv(bs_comp,"bs_comp.csv")
bs_comp


############################################################
#combine plots for Figure 3
library(cowplot)

#combine
f2a = R2.plot
f2b = plot_grid(act.p1,loc.p1, align = "h", nrow = 1)
f2b = plot_grid(f2b, bs.p1, nrow = 2)
figure2 = plot_grid(f2a, f2b, ncol = 1, rel_heights = c(0.35,0.65))
ggsave("figure2_primateivso.png", figure2, width = 11, height = 9)

############################################################
#habitat effects
datal2 = lapply(datal, FUN = function(x) {
  x = x[!is.na(x$Habitat_cat),]
  return(x)})

SO_m = bf(SO_counts | trials(SO_tot) ~
            1 +  mo(Habitat_heterogenity) + Habitat_cat + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 + mo(Habitat_heterogenity) + Habitat_cat + effort_wf +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.hab = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
           family = c(multinomial, binomial), data = datal2, data2 = phylol,
           prior = gen_prior,
           backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
           warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
           control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.hab, "m_hab.RDS")
m.hab = readRDS("m_hab.RDS")
summary(m.hab, robust = TRUE)
fixef(m.hab, robust = TRUE)

############################################################
#fruit diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~
            1 +  fruitprop + effort_wf +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
            1 +  fruitprop + effort_wf +
            (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.fruit = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
              family = c(multinomial, binomial), data = datal, data2 = phylol,
              prior = gen_prior,
              backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
              warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
              control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.fruit, "m_fruit.RDS")
m.fruit = readRDS("m_fruit.RDS")
summary(m.fruit, robust = TRUE)
fixef(m.fruit, robust = TRUE)

############################################################
#folivory diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~
            1 +  folivprop + effort +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 +  folivprop + effort +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.foliv = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
              family = c(multinomial, binomial), data = datal, data2 = phylol,
              prior = gen_prior,
              backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
              warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
              control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.foliv, "m_foliv.RDS")
m.foliv = readRDS("m_foliv.RDS")
summary(m.foliv, robust = TRUE)
fixef(m.foliv, robust = TRUE)


############################################################
#seed diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~
              1 +  seedprop + effort +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 +  seedprop + effort +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.seed = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
                      family = c(multinomial, binomial), data = datal, data2 = phylol,
                      prior = gen_prior,
                      backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
                      warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
                      control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.seed, "m_seed.RDS")
m.seed = readRDS("m_seed.RDS")
summary(m.seed, robust = TRUE)
fixef(m.seed, robust = TRUE)

############################################################
#animal diet proportion

SO_m = bf(SO_counts | trials(SO_tot) ~
            1 +  animalprop + effort +
            (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
              1 +  animalprop + effort +
              (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.animal = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
               family = c(multinomial, binomial), data = datal, data2 = phylol,
               prior = gen_prior,
               backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
               warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
               control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.animal, "m_animal.RDS")
m.animal = readRDS("m_animal.RDS")
summary(m.animal, robust = TRUE)
fixef(m.animal, robust = TRUE)

############################################################
#foraging style
SO_m = bf(SO_counts | trials(SO_tot) ~
    1 + foraging_style + effort +
    (1|superfamily) + (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|obs))
IVSO_m = bf(IVSOint | trials(SO_tot) ~
    1 + foraging_style + effort +
    (1|gr(phylo, cov = A)) + (1|Genus_species) + (1|superfamily) + (1|obs))

m.forg = brm_multiple(formula = SO_m + IVSO_m + set_rescor(FALSE),
             family = c(multinomial, binomial), data = datal, data2 = phylol,
             prior = gen_prior,
             backend="cmdstanr", stan_model_args=list(stanc_options = list("O1")),
             warmup = n_warm, iter=n_iter, chains = n_chains, init = 0,
             control=list(adapt_delta=0.95, max_treedepth=10))

saveRDS(m.forg, "m_forg.RDS")
m.forg = readRDS("m_forg.RDS")
summary(m.forg, robust = TRUE)
fixef(m.forg, robust = TRUE)

############################################################
#please email jordanscottmartin@gmail.com if you have any
#questions or comments



