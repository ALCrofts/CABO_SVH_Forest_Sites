# Author: Anna L. Crofts, email: croftsanna@gmail.com

# Title: Code for the manuscript, "Linking aerial hyperspectral data to canopy
#        tree biodiversity: an examination of the spectral variation hypothesis."

# Notes: 
#   Section 0: Load packages and data, etc.
#   Section 1: Degree of correspondence analyses 
#       Section 1.1: Composition - procrustes and pearson's correlations
#       Section 2.2: Diversity   - pearson's correlations
#   Section 2: Environmental drivers of community properties
#       Section 2.1: Composition - RDAs
#       Section 2.2: Diversity   - GAMs
#   Section 3: Extra Figures
#       Section 3.1: Study Site Figure
#       Section 3.2: Spurious Northness Effects
#       Section 3.3: Wavelengths Contributing to Spectral Diversity

# 0.0_Load Packages and data, etc. --------------------------------------

# Load Packages
library(tidyverse)
library(vegan)
library(hillR)
library(ggplot2)
library(ggtext)
library(stats) 
library(FD)
library(hillR)
library(corrplot)
library(ggpubr)
library(mgcv)
library(mgcViz)
library(MuMIn)
library(rstatix)

# Load Data
TaxComp <- read.csv("Data/VisibleCommunity_Corrected.csv")
FunTraits <- read.csv("Data/Trait_Species.csv")
CasiComp <- read.csv("Data/3_CASI_Processed_Spectra.csv")
SasiComp <- read.csv("Data/3_SASI_Processed_Spectra.csv")
Envr <- read.csv("Data/PlotEnvironment.csv")

# Define plots, species, and sites
plot_field_id <- TaxComp$plot_field_id

SpeciesList <- c("Abies.balsamea..Linnaeus..Miller", "Acer.pensylvanicum.Linnaeus", "Acer.rubrum.Linnaeus",                  
                 "Acer.saccharum.Marshall",          "Acer.spicatum.Lamarck",       "Alnus.incana..Linnaeus..Moench",        
                 "Betula.alleghaniensis.Britton",    "Betula.papyrifera.Marshall",  "Betula.populifolia.Marshall",           
                 "Carya.cordiformis..Wangenheim..K..Koch", "Fagus.grandifolia.Ehrhart",  "Fraxinus.americana.Linnaeus",           
                 "Fraxinus.nigra.Marshall",          "Larix.laricina..Du.Roi..K..Koch",  "Ostrya.virginiana..Miller..K..Koch",    
                 "Picea.abies..Linnaeus..H..Karsten",      "Picea.glauca..Moench..Voss", "Picea.rubens.Sargent",                  
                 "Pinus.strobus.Linnaeus",           "Populus.balsamifera.Linnaeus",     "Populus.grandidentata.Michaux",         
                 "Populus.tremuloides.Michaux",      "Prunus.pensylvanica.Linnaeus.f.",  "Prunus.serotina.Ehrhart.var..serotina", 
                 "Quercus.rubra.Linnaeus",           "Sorbus.decora..Sargent..C.K..Schneider", "Thuja.occidentalis.Linnaeus",           
                 "Tilia.americana.Linnaeus",         "Tsuga.canadensis..Linnaeus..Carrière",   "Ulmus.rubra.Muhlenberg") 

SiteList <- Envr %>%
  dplyr::select(site, field_plot_id) %>%
  rename(plot_field_id = field_plot_id)

# Calculate percent coniferous cover per plot
PercConifer <- TaxComp %>%
  group_by(plot_field_id) %>%
  summarise(PercCon = sum(c_across(cols = c('Abies.balsamea..Linnaeus..Miller',
                                            'Larix.laricina..Du.Roi..K..Koch',
                                            'Picea.abies..Linnaeus..H..Karsten',
                                            'Picea.glauca..Moench..Voss',
                                            'Picea.rubens.Sargent',
                                            'Pinus.strobus.Linnaeus',
                                            'Thuja.occidentalis.Linnaeus',
                                            'Tsuga.canadensis..Linnaeus..Carrière')))) %>%
  left_join(SiteList)

# 1.0_Degree of Correspondence Analyses  ---------------------------------------------------------
#   The spectral variation hypothesis predicts that variation in plant reflectance is related to 
#   variation in plant taxonomic and functional identity (Palmer et al. 2002, Ustin & Gamon 2010).
#   Spectral and field based community properities, both composition and diversity, should be related.

# 1.1__Composition ---------------------------------------------------------
# Examine the degree of correspondence of composition using two methods:
#   i)  Multivariate correlation - ie. Procrustes analysis
#   ii) Univariate correlation   - ie. Pearson's correlation of principal component axes

# 1.1.1__Procrustes Analyses ---------------------------------------------
# Procrustes analysis is a canonical ordination method for comparing two data matrices about the same objects (i.e., plots). Unlike Mantel test, 
# it is appropriate to use on raw data (vs. dissimilarity matrices). The procrustes test finds a compromise ordination for two raw data matrices (by dilating, rotating, and translating),
# where the sums of squared deviations are minimized. The symmetric procrustes statistic m12^2 is a measure of similarity between the two data matrices and can be tested by the 
# Procrustean randomization test.
#     - We want to compare the PC axes that explain most of the variation in compostition (>= 90% , ie, 'total' composition)
#     - The higher value of m12^2, the weaker the relationship between the two data matrices
#     - Procrustean r = (1-SS)^1/2


# 1.1.1.1__Explore Tax, Fun, and Spectral Comp ----------------------------

# 1.1.1.1.1__Taxonomic Composition ----------------------------------------

# Check the structure of species abundance data frame 
sum(TaxComp == 0) # 1633 zeros
sum(TaxComp == 0)/(nrow(TaxComp)*ncol(TaxComp)) # ~78% zeros

# There is quite a lot of zeros, so apply Hellinger Transformation 
#          -this transformation gives low weights to variables with low counts and many zeros
#          -square root of the relative abundance 
TaxComp.hellinger <- sqrt(TaxComp[,3:32])  # TaxComp is already characterized using relative abundances, so just need to square root transform.

# Ordinate the data and examine explained variance
#           -First two PCs explain 61.973% of the variance
#           -PC1 48.58%, PC2 13.394%, PC3 9.174% 
#           -90% variance explained by the PC1-PC8
pca.Tax <- rda(TaxComp.hellinger)
summary(pca.Tax) 

# Explore ordination
plot(pca.Tax, scaling = 1, display = "sites", type = "text", main = "PCA for Taxonomic Composition")

#   Extract species scores and plot them along PCs
pca.species.scores <- data.frame(scores(pca.Tax, display = 'species')) %>%
  rownames_to_column('Species')

PC1.Tax <- ggplot(pca.species.scores, aes(x = reorder(Species, PC1), y = PC1)) +
  geom_bar(stat = 'identity') +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Tax # Splits low and high elevation plots 

PC2.Tax <- ggplot(pca.species.scores, aes(x = reorder(Species, PC2), y = PC2)) +
  geom_bar(stat ='identity') +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Tax # Splits (kinda) Mt St Bruno and Mt Megantic


# 1.1.1.1.2__Functional Composition ----------------------------------------
# Calculate CWM of traits for each plot
FunTraits <- FunTraits %>%
  dplyr::select(-c("X", "scientific_name", "trait_leaf_mass_per_area_g_m2", "trait_actual_leaf_dry_matter_content_perc",
            "trait_leaf_water_content_mg_g", "trait_ndf_perc", "trait_adf_perc", "trait_adl_perc", 
            "trait_carot_mg_g_disk_mass", "trait_chlb_mg_g_disk_mass",  "trait_chla_mg_g_disk_mass",
            "trait_chlb_mg_l", "trait_carot_mg_l")) %>%
  add_column(SpeciesList, .before = 'trait_specific_leaf_area_m2_kg') 
  

FunComp <- TaxComp %>%
  dplyr::select(-'X') %>%
  gather(species, abundance, Abies.balsamea..Linnaeus..Miller:Ulmus.rubra.Muhlenberg) %>%
  group_split(plot_field_id) %>%
  map(~.x %>% left_join(FunTraits, by=c('species' = 'SpeciesList')), data=.x) %>%
  map_df(~.x %>% summarise_at(vars("trait_specific_leaf_area_m2_kg":"trait_c_perc"), ~weighted.mean(.x, abundance))) %>%
  add_column(plot_field_id, .before = ('trait_specific_leaf_area_m2_kg'))

# Standardize CWMs, so traits are scaled to zero mean and unit variance
FunComp.stand <- decostand(FunComp[,2:16], method = "standardize")

# Ordinate the stanadardized CWMs and examine explained variance
#           - First two PCs explain 78.783% of the variance
#           - PC1 57.76%, PC2 21.02%, PC3 14.44% 
#           - >99% variance explained by the PC1-PC8
pca.Fun <- rda(FunComp.stand)
summary(pca.Fun) 

# Explore ordination
plot(pca.Fun, scaling = 1, display = "sites", type = "text", main = "PCA for Functional Composition")

#   Extract trait scores and plot them along PCs
pca.trait.scores <- data.frame(scores(pca.Fun, display = 'species')) %>%
  rownames_to_column('Traits')

PC1.Fun <- ggplot(pca.trait.scores, aes(x = reorder(Traits, PC1), y = PC1)) +
  geom_bar(stat ='identity') +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Fun # Splits conifer and deciduous chemically - cellulose perc, lignin perc, recalcitrants perc, c perc,  chla:chlb, and equivalent H2O thickness  

PC2.Fun <- ggplot(pca.trait.scores, aes(x = reorder(Traits, PC2), y = PC2)) +
  geom_bar(stat='identity') +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Fun #Splits (kinda) conifer and deciduous structurally


# 1.1.1.1.3__VNIR-Spectral Composition ----------------------------------------
# Spectral composition is defined as the mean normalized reflectance across all bands (ie. wavelengths) per plot
CasiComp.mean <- CasiComp %>% 
  select(c(6,20:248)) %>%  # Spectral reflectance is in columns 20:248
  group_by(plot_field_id) %>%
  summarise_all(mean)
  
# Explore ordination
#           - First two PCs explain 97.67% of the variance
#           - PC1 84.63%, PC2 13.04%, PC3 0.96% 
#           - >99% variance explained by the PC1-PC8
pca.Casi <- rda(CasiComp.mean[,-1])
summary(pca.Casi) 

#   Extract wavelength scores and plot them along PCs
Casi.band.scores <- data.frame(pca.Casi$CA$v) %>%
  select(c(1:8)) %>%
  rownames_to_column(var="Bands") %>%
  mutate(Wavelength = gsub('BD_Band', "", Bands),
         Wavelength = gsub('nm', '', Wavelength),
         Wavelength = as.numeric(Wavelength))

PC1.Casi <- ggplot(Casi.band.scores, aes(x = Wavelength, y = PC1)) +
  geom_line() +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Casi 

PC2.Casi <- ggplot(Casi.band.scores, aes(x = Wavelength, y = PC2)) +
  geom_line() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Casi 


# 1.1.1.1.3__SWIR-Spectral Composition ----------------------------------------
SasiComp.mean <- SasiComp %>% 
  select(c(6,22:99)) %>%  # Spectral reflectance is in columns 20:248
  group_by(plot_field_id) %>%
  summarise_all(mean)

# Explore ordination
#           - First two PCs explain 86.93% of the variance
#           - PC1 58.68%, PC2 28.25%, PC3 7.23% 
#           - >99% variance explained by the PC1-PC8
pca.Sasi<-rda(SasiComp.mean[,-1])
summary(pca.Sasi) 

#   Extract wavelength scores and plot them along PCs
Sasi.band.scores <- data.frame(pca.Sasi$CA$v) %>%
  select(c(1:8)) %>%
  rownames_to_column(var="Bands") %>%
  mutate(Wavelength = gsub('BD_Band', "", Bands),
         Wavelength = gsub('nm', '', Wavelength),
         Wavelength = as.numeric(Wavelength))

PC1.Sasi <- ggplot(Sasi.band.scores, aes(x = Wavelength, y = PC1)) +
  geom_line () +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Sasi 

PC2.Sasi <- ggplot(Sasi.band.scores, aes(x = Wavelength, y = PC2)) +
  geom_line() +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Sasi


# 1.1.1.2_Run Procrustes Analyses -----------------------------------------

# PCs 1-8 contains greater than 90% of variation in composition across all dimensions 
#   Create dataframes of site scores across PCs 1-8

TaxPCA.scores <- data.frame(pca.Tax$CA$u) %>% 
  cbind(plot_field_id) %>%
  select(c('plot_field_id','PC1','PC2','PC3','PC4', 'PC5','PC6', 'PC7', 'PC8'))

FunPCA.scores <- data.frame(pca.Fun$CA$u) %>% 
  cbind(plot_field_id) %>%
  select(c('plot_field_id','PC1','PC2','PC3','PC4', 'PC5','PC6', 'PC7', 'PC8'))

CasiPCA.scores <- data.frame(pca.Casi$CA$u) %>% 
  cbind(plot_field_id) %>%
  select(c('plot_field_id','PC1','PC2','PC3','PC4', 'PC5','PC6', 'PC7', 'PC8'))

SasiPCA.scores <- data.frame(pca.Sasi$CA$u) %>% 
  cbind(plot_field_id) %>%
  select(c('plot_field_id','PC1','PC2','PC3','PC4', 'PC5','PC6', 'PC7', 'PC8'))

# Create and save dataframe with PC site scores for all dimensions

TaxPCA.scores2 <- TaxPCA.scores %>%
  mutate(Dimension = 'Taxonomic')

FunPCA.scores2 <- FunPCA.scores %>%
  mutate(Dimension = 'Functional')

CasiPCA.scores2 <- CasiPCA.scores %>%
  mutate(Dimension = 'Casi')

SasiPCA.scores2 <- SasiPCA.scores %>%
  mutate(Dimension = 'Sasi')

PCA_outputs <- rbind(TaxPCA.scores2, FunPCA.scores2, CasiPCA.scores2, SasiPCA.scores2)
write.csv(PCA_outputs, "Outputs/Data/PCA_outputs.csv")

rm(TaxPCA.scores2, FunPCA.scores2, CasiPCA.scores2, SasiPCA.scores2)


# Run procrustes analyses and procrustean permutation test (n = 999)
Tax.Casi <- procrustes(X = TaxPCA.scores[,2:9], Y = CasiPCA.scores[,2:9], symmetric = T)
Tax.Casi.out <- protest(X = TaxPCA.scores[,2:9], Y = CasiPCA.scores[,2:9], symmetric = T, scores = "sites", permutations = 999)

Tax.Sasi <- procrustes(X = TaxPCA.scores[,2:9], Y = SasiPCA.scores[,2:9], symmetric = T)
Tax.Sasi.out <- protest(X = TaxPCA.scores[,3:9], Y = SasiPCA.scores[,3:9], symmetric = T, scores = "sites", permutations = 999)

Fun.Casi <- procrustes(X = FunPCA.scores[,2:9], Y = CasiPCA.scores[,2:9], symmetric = T)
Fun.Casi.out <- protest(X = FunPCA.scores[,2:9], Y = CasiPCA.scores[,2:9], symmetric = T, scores = "sites", permutations = 999)

Fun.Sasi <- procrustes(X = FunPCA.scores[,2:9], Y = SasiPCA.scores[,2:9], symmetric = T)
Fun.Sasi.out <- protest(X = FunPCA.scores[,2:9], Y = SasiPCA.scores[,2:9], symmetric = T, scores = "sites", permutations = 999)

# Create a dataframe detailing procustes output
Procrustes <- data.frame(Field = c('Taxonomic', 'Taxonomic', 'Functional', 'Functional'),
                         Spectral = c('Casi', 'Sasi', 'Casi', 'Sasi'),
                         r =c(Tax.Casi.out$t0, Tax.Sasi.out$t0, Fun.Casi.out$t0, Fun.Sasi.out$t0),
                         m12.squared =c(Tax.Casi.out[["ss"]], Tax.Sasi.out[["ss"]], Fun.Casi.out[["ss"]], Fun.Sasi.out[["ss"]]),
                         p = c(Tax.Casi.out[["signif"]], Tax.Sasi.out[["signif"]], Fun.Casi.out[["signif"]], Fun.Sasi.out[["signif"]]))

write.csv(Procrustes, "Outputs/Statistics/Procrustes.csv")

# Plot fit of procrustes analyses  
ProcrustesFit<-ggplot(Procrustes, aes(x = Spectral, y = r, fill = factor(Field, levels =c("Taxonomic", "Functional")))) +
  geom_bar(stat = 'identity', position = position_dodge(), colour = 'black') +
  scale_fill_manual(values = c("Grey40", "Grey80")) +
  #geom_text(aes(label = r), position=position_dodge(width=0.9), vjust=-0.25, size=4)+
  scale_y_continuous(expand=c(0,0),breaks=seq(0, 1, by=0.25), limits=c(0,1))+ 
  scale_x_discrete(labels=c('VNIR', 'SWIR'))+
  theme_classic() +
  ylab("Procrustean r")+
  xlab("")+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(face='bold', size=12, colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(face='bold', size=10),
        axis.title.y =element_text(face='bold', size=12),
        legend.text = element_text(face='bold', size=12),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.95)) 
ProcrustesFit
ggsave("Outputs/Figures/ProcrustesFit.jpeg", width=6, height=4, dpi=300)

# # Run procrustes again but with just the 1st three PCs (vs. PCs explaining >90% variance)
# Tax.Casi3 <- procrustes(X = TaxPCA.scores[,2:4], Y = CasiPCA.scores[,2:4], symmetric = T)
# Tax.Casi.out3 <- protest(X = TaxPCA.scores[,2:4], Y = CasiPCA.scores[,2:4], symmetric = T, scores = "sites", permutations = 999)
# 
# Tax.Sasi3 <- procrustes(X = TaxPCA.scores[,2:4], Y = SasiPCA.scores[,2:4], symmetric = T)
# Tax.Sasi.out3 <- protest(X = TaxPCA.scores[,2:4], Y = SasiPCA.scores[,2:4], symmetric = T, scores = "sites", permutations = 999)
# 
# Fun.Casi3 <- procrustes(X = FunPCA.scores[,2:4], Y = CasiPCA.scores[,2:4], symmetric = T)
# Fun.Casi.out3 <- protest(X = FunPCA.scores[,2:4], Y = CasiPCA.scores[,2:4], symmetric = T, scores = "sites", permutations = 999)
# 
# Fun.Sasi3 <- procrustes(X = FunPCA.scores[,2:4], Y = SasiPCA.scores[,2:4], symmetric = T)
# Fun.Sasi.out3 <- protest(X = FunPCA.scores[,2:4], Y = SasiPCA.scores[,2:4], symmetric = T, scores = "sites", permutations = 999)
# 
# Procrustes_PC3 <- data.frame(Field = c('Taxonomic', 'Taxonomic', 'Functional', 'Functional'),
#                          Remote = c('Casi', 'Sasi', 'Casi', 'Sasi'),
#                          r =c(Tax.Casi.out3$t0, Tax.Sasi.out3$t0, Fun.Casi.out3$t0, Fun.Sasi.out3$t0),
#                          m12.squared =c(Tax.Casi.out3[["ss"]], Tax.Sasi.out3[["ss"]], Fun.Casi.out3[["ss"]], Fun.Sasi.out3[["ss"]]),
#                          p = c(Tax.Casi.out3[["signif"]], Tax.Sasi.out3[["signif"]], Fun.Casi.out3[["signif"]], Fun.Sasi.out3[["signif"]]))

rm(Tax.Casi, Tax.Casi.out, Tax.Sasi, Tax.Sasi.out, Fun.Casi, Fun.Casi.out, Fun.Sasi, Fun.Sasi.out, ProcrustesFit, Procrustes)

# 1.1.2_Pearson's Correlations ----------------------------------------------------
#   Examine if the individual PCs are associated with each other.
#   Visualize species, trait, wavelenght scores along PC axes


# 1.1.2.1__Run PC Correlations --------------------------------------------

# PCA_outputs <- read.csv("Outputs/Data/PCA_outputs.csv")

PCA_outputs2 <- PCA_outputs %>%
  #dplyr::select(-c('X')) %>%
  pivot_wider(names_from = Dimension, values_from = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8))

# Pearson's correlations of the first 3 PC axes
PCACor.Total <- cor(PCA_outputs2[,2:13])
colnames(PCACor.Total) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                            "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                            "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")
rownames(PCACor.Total) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                            "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                            "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")

PCACor.Test = cor.mtest(PCA_outputs2[,2:13], conf.level = 0.95)
colnames(PCACor.Test[['p']]) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                                  "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                                  "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")
rownames(PCACor.Test[['p']]) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                                  "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                                  "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")

# # Create a dataframe with correlation coefficients for the first 3 PCs across dimensions
PCACor.Total2 <- pull_lower_triangle(PCACor.Total) %>% 
  gather(Var2, r, 2:13) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) %>%
  mutate(r = round(r, 2)) 

PCACor.Test2 <- pull_lower_triangle(PCACor.Test[["p"]]) %>% 
  gather(Var2, p, 2:13) %>%
  rename(Var1 = rowname) %>%
  mutate(p = as.numeric(p)) %>%
  filter(!is.na(p)) %>%
  mutate(p = round(p, 3)) 

PCACor.Total3 <- left_join(PCACor.Total2, PCACor.Test2) %>%
  filter(Var1 == "VNIR_PC1" | Var1 == "VNIR_PC2" | Var1 == "VNIR_PC3" |
         Var1 == "SWIR_PC1" | Var1 == "SWIR_PC2" | Var1 == "SWIR_PC3" ) %>%
  filter(Var2 == "Tax_PC1" | Var2 == "Tax_PC2" | Var2 == "Tax_PC3" |
         Var2 == "Fun_PC1" | Var2 == "Fun_PC2" | Var2 == "Fun_PC3" ) %>%
  rename(Spectra = Var1, Field = Var2)
write.csv(PCACor.Total3, "Outputs/Statistics/PCA_correlations.csv")

# Plot the Correlation matrix 
PCACor.Total4 <- left_join(PCACor.Total2, PCACor.Test2) %>%
  mutate(r2 = ifelse(p <= 0.05, r, NA))  # Remove r values for non-signif correlations -- FOR PLOTTING PURPOSES!!

PCACorTotalMatrix <- ggplot(data = PCACor.Total4, aes(Var2, Var1, fill = r))+
  scale_x_discrete(limits = c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                              "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                              "Tax_PC3", "Fun_PC3", "VNIR_PC3")) +
  scale_y_discrete(limits = c("SWIR_PC3","VNIR_PC3","Fun_PC3", "Tax_PC3",
                              "SWIR_PC2","VNIR_PC2","Fun_PC2", "Tax_PC2",
                              "SWIR_PC1","VNIR_PC1","Fun_PC1")) +
  geom_tile(color = "white")+
  labs(title = "") +
  geom_text(aes(Var2, Var1, label = r), color = "black", size = 4) +
  geom_text(aes(Var2, Var1, label = r2), na.rm = T, colour = 'black', fontface = 'bold', size = 4) +
  scale_fill_distiller(palette = 'RdBu', direction = 1, limit = c(-1,1), name="Pearson's Correlation",
                       guide = guide_colourbar(title.hjust = 0.5, title.position = "top")) +
  theme_void()+ 
  theme(axis.text.x = element_text(face='bold', colour = 'black', angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(face='bold', size=12, colour = 'black'),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_text(face='bold', size = 12),
        #legend.text.align = 0.5,
        legend.position = 'bottom',
        #legend.box.background = element_rect(colour = "black"),
        #legend.title.align=0.5,
        legend.key.width = unit(2, "cm")) +
  coord_fixed() 
PCACorTotalMatrix

ggsave(plot = PCACorTotalMatrix, filename = 'Outputs/Figures/PCACor_Total.jpeg', width=6, height=6, dpi=300)

rm(PCACor.Total, PCACor.Total2, PCACor.Total3, PCACor.Total4, PCACor.Test, PCACor.Test2, PCACorTotalMatrix)

# 1.1.2.2__Understand PC axes -------------------------------------
#   - How do the primary principal components relate to the temperate-to-boreal gradient? 
#   - How are species, traits, and wavelengths distributed along the primary principal components?

# Pearson's correlations of primary PC axes and percent conifer content
PCA_outputs3 <- PCA_outputs2 %>%
  left_join(PercConifer)

PCACor.PC1 <- cor(PCA_outputs3[,c(2:5,34)])
colnames(PCACor.PC1) <- c("Taxonomic", "Functional", "VNIR-Spectral", "SWIR-Spectral", "Perc. Coniferous")
rownames(PCACor.PC1) <- c("Taxonomic", "Functional", "VNIR-Spectral", "SWIR-Spectral", "Perc. Coniferous")

PCACor.PC1_2 <- pull_lower_triangle(PCACor.PC1) %>% 
  gather(Var2, r, 2:6) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) %>%
  mutate(r = round(r, 2)) 

PC1CorMatrix <- ggplot(data = PCACor.PC1_2, aes(Var2, Var1, fill = r))+
  scale_x_discrete(limits = c('Taxonomic','Functional','VNIR-Spectral','SWIR-Spectral')) +
  scale_y_discrete(limits = c('Perc. Coniferous', 'SWIR-Spectral', 'VNIR-Spectral', 'Functional')) +
  geom_tile(color = "white")+
  labs(title = "") +
  geom_text(aes(Var2, Var1, label = r), color = "black", size = 5) +
  scale_fill_distiller(palette = 'RdBu', direction = 1, limit = c(-1,1), name="Pearson's Correlation",
                       guide = guide_colourbar(title.hjust = 0.5, title.position = "top")) +
  theme_minimal()+ 
  theme(axis.text.x = element_text(face='bold', colour = 'black', angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
      axis.ticks.x=element_blank(),
      axis.text.y = element_text(face = 'bold', size=12, colour = 'black'),
      axis.title.y=element_blank(),
      axis.title.x = element_blank(),
      legend.title=element_text(face='bold', size = 12),
      #legend.text.align = 0.5,
      legend.position = 'bottom',
      #legend.box.background = element_rect(colour = "black"),
      #legend.title.align=0.5,
      legend.key.width = unit(1, "cm")
      ) +
  coord_fixed() 
PC1CorMatrix

# Plot highly correlated PC axes: PC 1 for Tax, Fun, CASI and PC 2 for SASI
# Filter species scores, so only 15 most abundant species retained 
TopSpecies <- TaxComp %>%
  select(-c(1:2)) %>%
  summarise_all(sum) %>%
  gather(key = 'Species', value = 'Abundance') %>%
  arrange(desc(Abundance)) %>%
  slice(1:15) 

pca.species.scores2 <- pca.species.scores %>%
  filter(Species %in% TopSpecies$Species) %>%
  mutate(Conifer = c('C', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'C', 'D', 'C', 'D', 'D', 'C')) %>%
  mutate(Names = c('Abies balsamea', 'Acer pensilvanicum', 'Acer rubrum', 'Acer saccharum', 'Betula alleganiensis', 'Betula papyrifera', 'Betula popufolia', 'Fagus grandifolia','Ostrya virginiana','Picea rubens', 'Populus tremuloides','Prunus strobus','Quercus rubra','Tilia americana','Tsuga canadensis'))

PC1.Tax <- ggplot(pca.species.scores2, aes(x=reorder(Names, PC1), y=PC1)) +
  geom_bar(stat='identity', colour = 'black', fill = 'grey60') +
  #scale_fill_manual(values = c('#004a00', '#8cd66c')) + 
  xlab("")+
#  ylab("Taxonomic PC1 \n (49%) ") + 
  ylab('PC1 Scores') +
  labs(title = 'Taxonomic', subtitle = 'Variance Explained = 49%') +
#  labs(title = 'Taxonomic') +
  scale_x_discrete(labels = c("ACSA", "FAGR", "QURU", "OSVI","TSCA","TIAM","PRSE","ACPE","BEPO","POTR","ACRU","BEAL","PIRU","BEPA","ABBA")) +
  scale_y_continuous(limits = c(-1.25, 1.25)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face='bold', size=12),
        plot.subtitle = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1, face='bold', size=10, color='black'),
        axis.text.y = element_text(face='bold', size=8, color='black'),
        axis.title.y = element_text(face='bold', size=10),
        axis.title.x = element_text(face='bold', size=12),
        legend.title = element_text(face='bold', size = 10)) 
PC1.Tax #Splits low and high elevation plots 

pca.trait.scores2 <- pca.trait.scores %>%
  mutate(Traits2 = c('SLA', 'LDMC','RWC', 'EWT', 'Soluable', 'Hemicellulose', 'Cellulose', 'Lignin', 'Recalcitrants', 'Chl a', 'Chl b', 'Carotenoids', 'Chl a : Chl b', 'N', 'C'))

PC1.Fun <- ggplot(pca.trait.scores2, aes(x=reorder(Traits2, PC1), y=PC1)) +
  geom_bar(stat='identity', colour = 'black', fill = 'grey60') +
  xlab("")+
#  ylab("Functional PC1 \n (58%)") +
  ylab("PC1 Scores") +
  scale_x_discrete(labels = c("Carot.", "SLA","Chl b", "Chl a", "RWC", "N", "Hemi.", "LDMC", "Sol.", "Cell.", "Chla:Chlb", "Lignin","Recalc.", "C", "EWT")) +
  labs(title = 'Functional', subtitle = 'Variance Explained = 58%') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face='bold', size=12),
        plot.subtitle = element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 1, face='bold', size=10, color='black'),
        axis.text.y = element_text(face='bold', size=8, color='black'),
        axis.title.y = element_text(face='bold', size=10),
        axis.title.x = element_text(face='bold', size=12),
        legend.title = element_text(face='bold', size = 10),
        legend.position = 'none') 
PC1.Fun #Splits conifer and deciduous chemically - cellulose perc, lignin perc, recalcitrants perc, c perc,  chla:chlb, and equivalent H2O thickness  

Casi.band.scores2 <- Casi.band.scores %>%
  mutate(Region = ifelse(Wavelength <= 500, 'blue',
                         ifelse(Wavelength > 500 & Wavelength <= 600, 'green',
                                ifelse(Wavelength > 600 & Wavelength <= 730, 'red',
                                       ifelse(Wavelength > 730 & Wavelength <=780, 'rededge', 'NIR'))))) %>%
  full_join(data.frame(Region = 'SWIR',
                       PC1 = 0,
                       Wavelength = 454))
PC1.Casi <- ggplot(Casi.band.scores2, aes(x=Wavelength, y = PC1)) +
  geom_area(fill = 'steelblue4') +
  geom_area(mapping = aes(x = ifelse(Wavelength > 590, Wavelength, 0)), fill = 'mediumseagreen') +
  geom_area(mapping = aes(x = ifelse(Wavelength > 700, Wavelength, 0)), fill = 'tomato3') +
  geom_area(aes(fill = Region)) +
  geom_line() +
  scale_fill_manual(values = c('steelblue4', 'mediumseagreen','tomato3','lightcoral',  'indianred3',  'rosybrown3'),
                    breaks = c('blue','green','red','rededge','NIR','SWIR'),
                    name = 'Spectral Region', 
                    labels = c('Blue','Green','Red', 'Red Edge', 'NIR', 'SWIR'),
                    guide = guide_legend(title.hjust = 0.5, title.position = "top", nrow = 1)) + 
  # geom_area(data = Casi.band.scores2 %>% filter(Region == 'green'), mapping = aes(x=Wavelength, y = PC1), fill ='green') +
  # geom_area(data = Casi.band.scores2 %>% filter(Region == 'red'), mapping = aes(x=Wavelength, y = PC1), fill ='red') +
  # geom_area(data = Casi.band.scores2 %>% filter(Region == 'rededge'), mapping = aes(x=Wavelength, y = PC1), fill ='darkred') +
  # geom_area(data = Casi.band.scores2 %>% filter(Region == 'NIR'), mapping = aes(x=Wavelength, y = PC1), fill ='purple') +
  annotate('rect', xmin = 899, xmax = 957, ymin = -0.21, ymax = 0.21, colour = 'grey90', fill = 'grey90') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab("Wavelength (nm)") + 
#  ylab("VNIR-Spectral PC1 \n (85%)") + 
  ylab('PC1 Scores') +
  labs(title = 'VNIR-Spectral', subtitle = 'Variance Explained = 84%') +
  scale_y_continuous(limits = c(-0.21, 0.21), expand = c(0,0)) +
  scale_x_continuous(limits = c(454,1059), breaks = seq(500,1000,100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(face='bold', size=12),
        plot.subtitle = element_text(size=10),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(face='bold', size=8, color='black'),
        axis.title.y = element_text(face='bold', size=10),
        axis.title.x = element_text(face='bold', size=10),
        legend.title = element_text(face='bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5)
PC1.Casi 

Sasi.band.scores2 <- Sasi.band.scores %>%
  mutate(Region = ifelse(Wavelength <= 1400, 'NIR', 'SWIR'))

PC2.Sasi <- ggplot(Sasi.band.scores2, aes(x=Wavelength, y=PC2, fill = Region)) +
  geom_area() +
  geom_line() +
  scale_fill_manual(values = c('indianred3', 'rosybrown3')) +
  annotate('rect', xmin=c(1333,1790), xmax=c(1466,1950), ymin=c(-0.33,-0.33), ymax=c(0.33,0.33), colour = 'grey90', fill = 'grey90') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab("Wavelength (nm)") +
#  ylab("SWIR-Spectral PC2 \n (26%)") +
  ylab("PC2 Scores") +
  labs(title = 'SWIR-Spectral', subtitle = 'Variance Explained = 28%') +
  scale_y_continuous(limits = c(-0.33, 0.33), breaks = seq(-0.3, 0.3, 0.15), expand = c(0,0)) +
  scale_x_continuous(limits = c(1002,2412), breaks = seq(1000,2400,200)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face='bold', size=12),
        plot.subtitle = element_text(size=10),
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(face='bold', size=8, color='black'),
        axis.title.y = element_text(face='bold', size=10),
        axis.title.x = element_text(face='bold', size=10),
        legend.title = element_text(face='bold', size = 10),
        legend.position = 'bottom') 
PC2.Sasi

windows(width = 7.08661, height = 6)
PCPlots <- ggarrange(PC1.Tax, PC1.Fun, PC1.Casi, PC2.Sasi,
          nrow = 2, ncol = 2,
          common.legend = T, legend = 'bottom',
          labels =  'AUTO')
PCPlots
ggsave(plot = PCPlots, filename = 'Outputs/Figures/PC_plots3.jpeg', width=7.08661, height=6, dpi=600)


rm(PCA_outputs, PCA_outputs2, PCA_outputs3, PCACor.PC1, PCACor.PC1_2, PC1CorMatrix,
   TaxComp.hellinger, TaxPCA.scores, pca.Tax, pca.species.scores,  pca.species.scores2, TopSpecies,
   FunComp.stand, FunPCA.scores,pca.Fun, pca.trait.scores, pca.trait.scores2,
   CasiComp.mean, CasiPCA.scores, pca.Casi, Casi.band.scores, Casi.band.scores2,
   SasiComp.mean, SasiPCA.scores, pca.Sasi, Sasi.band.scores, Sasi.band.scores2,
   PC1.Tax, PC1.Fun, PC1.Casi, PC1.Sasi, PC2.Tax, PC2.Fun, PC2.Casi, PC2.Sasi, PCPlots)

# 1.2__Diversity  --------------------------------------------------------------
#   - First, calculate diversity indices for each biodiversity dimension
#   - Then, examine the degree of association between spectral and taxonomic/functinal indices
#       * 1st via linear correlation for all indices 
#       * 2nd via non-linear association for subset (?) - best correlated indices


# 1.2.1__Calculate Diversity Indices --------------------------------------


# 1.2.1.1__Taxonomic Diversity Indices  -------------------------------------------------------------------------
# Calculate hill diversity for each field plot:
#   q0 = species richness, the number of species per plot
#   q1 = Shannon diversity, the uncertainty in predicting the identity of a new indiv. given species richness and abundance
#   q2 = Simpson diversity, the prob. that two indiv. drawn at random will be different species
#   J  = Pielou's eveness, how close in abundance each species is

TaxDiv <- TaxComp %>%
  summarise(plot_field_id = plot_field_id,
            q0 = hill_taxa(TaxComp[,3:32], q = 0, MARGIN = 1), 
            q1 = hill_taxa(TaxComp[,3:32], q = 1, MARGIN = 1), 
            q2 = hill_taxa(TaxComp[,3:32], q = 2, MARGIN = 1)) %>%
  mutate(J =  log(q1)/log(q0)) %>%
  mutate(q1 = ifelse(plot_field_id == 'Plot13_Maples', NA, q1),
         q2 = ifelse(plot_field_id == 'Plot13_Maples', NA, q2),
         J  = ifelse(plot_field_id == 'Plot13_Maples', NA, J))

# 1.2.1.2__Functional Diversity Indices ----------------------------------------
# Use FD package to calculate functional diversity indices for each plot:
#   FRic = Functional Richness, convex hull volume of functional trait space
#   FEve = Functional Evenness, quantifies the regularity of interspecies distances and the homogeneity of species abundances (btwn 0-1)
#   FDiv = Functional Divergence, the deviation of species from the centre of gravity functional trait space
#   FDis = Functional Dispersion, the mean distance of a species from the centroid of all species in functional trait space
#   RaoQ = Rao's Quadratic Entropy, the sum of distances of btwn species pairs wighted by their relative abundance 


dbFD <- dbFD(FunTraits %>% column_to_rownames(var='SpeciesList'), TaxComp[,3:32], w.abun = T, stand.x = T, calc.FRic = T, calc.FDiv = T)

FunDiv <- data.frame(plot_field_id = plot_field_id,
                     FRic = dbFD[["FRic"]],
                     FEve = dbFD[["FEve"]],
                     FDis = dbFD[['FDis']],
                     FDiv = dbFD[['FDiv']],
                     RaoQ = dbFD[['RaoQ']])

# 1.2.1.3__VNIR-Spectral Diversity Indices ----------------------------------------------
# Calculate three spectral diversity metrics for each plot: 
#   CV  = coefficient of variation, the average coefficient of variation for all wavelengths of the data within each plot. Larger CV corresponds to higher SD.
#   CHV = convex hull volume, the volume of a polygon of pixels/points forming a convex hull of the first three PC. Larger CHV corresponds with higher SD.
#   SV = spectral variance, the total sums of squared deviation for a community standardized by number of spectral points in the community (Laliberte et al. 2020). Higher SDa (aka neighboring pixels are spectral different from one another) corresponds with greater diversity. 

bootstrap <- function(df, n, p){  # n = the number of random sampling events, p = number of spectral points to randomly sample
  x <- replicate(n, sample_n(df, p, replace = F), F) # replace = F means that we are sampling without replacement
}

calc.CV <- function(df, grouping_var, .x){
  df %>% group_by({{grouping_var}}) %>%
    summarise_all(~sd(.x)/mean(.x)) %>%
    rowwise({{grouping_var}}) %>%
    summarise(CV = sum(c_across(cols = everything()), na.rm = T) / (ncol(.) - sum(is.na(c_across(everything())))))
}

calc.SV <- function(df, grouping_var, .x){
  
  points <- df %>%
    group_by({{grouping_var}}) %>%
    summarise(points = n())
  
  SV.df  <- df %>% group_by({{grouping_var}}) %>%
    summarise_all(~sum((.x - mean(.x))^2)) %>%
    rowwise({{grouping_var}}) %>%
    summarise(SS = sum(c_across(cols = everything()))) %>%
    left_join(points) %>%
    summarise(SV = SS / (points - 1))
  return(SV.df)
}

#  Create list of resampled dataframes 
Casi.boot <- CasiComp %>%
  select(c(6,20:248)) %>%
  group_by(plot_field_id) %>%
  bootstrap(999, 110) # p = 110 b/c this is the minimum number of spectral points observed in a plot

# Calculate bootstraped CV for each plot 
Casi.CV <- Casi.boot %>%
  map(~calc.CV(.x, plot_field_id, c(2:230))) %>%
  do.call(rbind, .) %>%
  group_by(plot_field_id) %>%
  summarise(CV = mean(CV))

# Calculate bootstraped CHV for each plot 
Casi.PCA.total <- CasiComp %>%
  select(c(20:248)) %>%
  rda(scale = F) 

Casi.CHV <- data.frame(Casi.PCA.total$CA$u) %>%
  select(c(1:3)) %>%
  cbind(plot_field_id = CasiComp$plot_field_id) %>%
  group_split(plot_field_id) %>%
  map(~bootstrap(.x, 999, 110) %>%
        map(~convhulln(.x[-4], option = 'FA')) %>%
        map_dbl('vol')) %>%
  map(~mean(.x)) %>%
  do.call(rbind, .) %>%
  cbind(plot_field_id = plot_field_id) %>%
  data.frame() %>%
  rename(CHV = V1)


# Calculate SV for each plot - SV is normalized per point number so no need to bootstrap
Casi.SV <- CasiComp %>%
  select(c(6,20:248)) %>%
  calc.SV(plot_field_id, c(2:230))

# Create CasiDiv dataframe
CasiDiv <- Casi.CV %>%
  left_join(Casi.CHV) %>%
  left_join(Casi.SV)
  
CasiDiv<- CasiDiv %>%
  rename(CV.Casi = CV,
         CHV.Casi = CHV,
         SV.Casi = SV)

#write.csv("Outputs/Data/CasiDiv.csv")
rm(Casi.boot, Casi.PCA.total)

# 1.2.1.4__SWIR-Spectral Diversity Indices ----------------------------------------------
# Create list of resampled dataframes 
Sasi.boot <- SasiComp %>%
  select(c(6,22:99)) %>%
  group_by(plot_field_id) %>%
  bootstrap(999, 95) # p = 95 b/c this is the minimum number of SASI spectral points observed in a plot

# Calculate bootstraped CV for each plot 
Sasi.CV <- Sasi.boot %>%
  map(~calc.CV(.x, plot_field_id, c(2:79))) %>%
  do.call(rbind, .) %>%
  group_by(plot_field_id) %>%
  summarise(CV = mean(CV))

# Calculate bootstraped CHV for each plot 
Sasi.PCA.total <- SasiComp %>%
  select(c(22:99)) %>%
  rda(scale = F) 

Sasi.CHV <- data.frame(Sasi.PCA.total$CA$u) %>%
  select(c(1:3)) %>%
  cbind(plot_field_id = SasiComp$plot_field_id) %>%
  group_split(plot_field_id) %>%
  map(~bootstrap(.x, 999, 95) %>%
        map(~convhulln(.x[-4], option = 'FA')) %>%
        map_dbl('vol')) %>%
  map(~mean(.x)) %>%
  do.call(rbind, .) %>%
  cbind(plot_field_id = plot_field_id) %>%
  data.frame() %>%
  rename(CHV = V1)

# Calculate SV for each plot - SV is normalized per point number so no need to bootstrap
Sasi.SV <- SasiComp %>%
  select(c(6,22:99)) %>%
  calc.SV(plot_field_id, c(2:79))

# Create SasiDiv dataframe
SasiDiv <- Sasi.CV %>%
  left_join(Sasi.CHV) %>%
  left_join(Sasi.SV)

SasiDiv<- SasiDiv %>%
  rename(CV.Sasi = CV,
         CHV.Sasi = CHV,
         SV.Sasi = SV)

rm(Sasi.boot, Sasi.PCA.total)


# 1.2.2__Diversity Correlations ----------------------------------------------
#   Examine the association of between spectral diversity indices and taxonomic / function diversity indices.
#     Pearson correlation measures the strength of the linear relationship between two variables

# First, create a dataframe containing all diversity indices
DivIndices <- TaxDiv %>%
  left_join(FunDiv) %>%
  left_join(CasiDiv) %>%
  left_join(SasiDiv) %>%
  mutate(CHV.Casi = as.numeric(CHV.Casi),
         CHV.Sasi = as.numeric(CHV.Sasi))

write.csv(DivIndices, file = "Outputs/Data/DiversityIndices.csv")

rm(Casi.CHV, Casi.CV, Casi.SV, CasiDiv,
   Sasi.CHV, Sasi.CV, Sasi.SV, SasiDiv,
   dbFD, FunDiv, TaxDiv)

# Do correlations and keep correlation statistics between spectral div indices and tax/fun div indices
#   - Pearson's correlation
#   - Complete.obs removes NA values by case wise deletion

# DivIndices <- read.csv('outputs/Data/DiversityIndices.csv') %>%
#   select(-c(1))

DivCor.total <- cor(DivIndices[,-1], use="complete.obs") #Do correlation, where we exclude Na values (ie. plots where div index could be calculated)
colnames(DivCor.total) <- c("q0", "q1", "q2", "J",
                            "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                            "CV.Casi", "CHV.Casi", "SD.Casi",
                            "CV.Sasi", "CHV.Sasi", "SD.Sasi")
rownames(DivCor.total) <- c("q0", "q1", "q2", "J",
                            "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                            "CV.Casi", "CHV.Casi", "SD.Casi",
                            "CV.Sasi", "CHV.Sasi", "SD.Sasi")

DivCor.Test = cor.mtest(DivIndices[,-1], use="complete.obs", conf.level = 0.95)
colnames(DivCor.Test[['p']]) <- c("q0", "q1", "q2", "J",
                                  "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                                  "CV.Casi", "CHV.Casi", "SD.Casi",
                                  "CV.Sasi", "CHV.Sasi", "SD.Sasi")
rownames(DivCor.Test[['p']]) <- c("q0", "q1", "q2", "J",
                                  "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                                  "CV.Casi", "CHV.Casi", "SD.Casi",
                                  "CV.Sasi", "CHV.Sasi", "SD.Sasi")

DivCor.Total2 <- pull_lower_triangle(DivCor.total) %>% 
  gather(Var2, r, 2:16) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) 

DivCor.Test2 <- pull_lower_triangle(DivCor.Test[["p"]]) %>% 
  gather(Var2, p, 2:16) %>%
  rename(Var1 = rowname) %>%
  mutate(p = as.numeric(p)) %>%
  filter(!is.na(p))  

DivCor <- left_join(DivCor.Total2, DivCor.Test2) %>%
  filter(Var1 == "CV.Casi" | Var1 == "CHV.Casi" |Var1 == "SD.Casi" |
           Var1 == "CV.Sasi" | Var1 == "CHV.Sasi" |Var1 == "SD.Sasi") %>%
  filter(Var2 == "q0" | Var2 == "q1" | Var2 == "q2" | Var2 == "J" |
           Var2 == "FRic" | Var2 == "FEve" | Var2 == "FDiv" | Var2 == "FDis" | Var2 == "RaoQ") %>%
  rename(Spectral = Var1, Field = Var2)

write.csv(DivCor, file = "Outputs/Statistics/DiversityCorrelations.csv")


# Summarize correlations by Sensor, Field Dimension, and Spectral Div metric and Plot
DivCor2 <- DivCor %>%
  separate(Spectral, c("Metric","Sensor"), sep = '[.]') %>%
  mutate(OnGround = ifelse(Field == "q0" |
                             Field == "q1" | 
                             Field == "q2" | 
                             Field == "J", "Tax", "Func")) 

DivCor2Sum <- DivCor2 %>%
  group_by(Sensor, OnGround, Metric) %>%
  mutate(mean = mean(r), se = sd(r)/sqrt(n()))

windows(width = 5, height = 4.5)
DivCorSum <- ggplot(DivCor2Sum, aes(x = factor(OnGround, level = c('Tax','Func')), y = mean, group = factor(Metric, level = c('CV','CHV',"SD")))) +
  geom_bar(stat = 'identity', position = position_dodge(), aes(fill = factor(Metric, level = c('CV', 'CHV', 'SD'))), colour = 'black') +
  geom_errorbar(aes(ymin = mean - se, ymax = mean +se), position = position_dodge(0.9), width = 0.2) +
  facet_wrap(~Sensor, labeller = labeller(Sensor = c("Casi" = "VNIR", "Sasi" = "SWIR"))) +
  scale_fill_grey(start = 0.8, end = 0.2, name = 'Spectral Metric', labels = c('CV','CHV','SV')) +
  scale_x_discrete("", labels = c("Taxonomic", 'Functional')) +
  scale_y_continuous("Mean Correlation with Spectral Diversity \n(Pearson's r)", expand = c(0,0), limits = c(-0,1), breaks = c(0,.25, .5,.75,1))+
  theme_classic() +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(axis.ticks.x=element_blank(),
        axis.text.x=element_text(face='bold', size=12, color='black'),
        axis.text.y = element_text(size=10, colour = 'black'),
        axis.title.y=element_text(face='bold', size=12),
        legend.title=element_text(face='bold', size = 12),
        legend.text.align = 0.5,
        legend.position = 'bottom',
        #legend.box.background = element_rect(colour = "black"),
        legend.title.align=0.5,
        strip.text = element_text(face = 'bold', size = 12))
DivCorSum

ggsave(plot = DivCorSum, filename = 'Outputs/Figures/DiversityCorrelations.jpeg', width=5, height=4.5, dpi=600)

# Create wide format
DivCor3 <- DivCor2 %>%
  mutate(r = round(r, 3)) %>%
  mutate(p = round(p, 3)) %>%
  pivot_wider(names_from = Metric, values_from = c("r","p")) 

write.csv(DivCor3, "Outputs/Statistics/DiversityCorrelations_wide.csv")

rm(DivCor, DivCor.Test, DivCor.Test2, DivCor.total, DivCor.Total2, DivCor2, DivCor3, DivCor2Sum, DivCorSum)

# 1.2.2.1__Plot Relationships w VNIR-Diversity --------------------------------------------------------------------
#   Further explore relationships between VNIR-spectral diversity and taxonomic/functional diversity

# Centre and standardize all spectral diversity indices and plot spectral ~ field 
DivInd.stand<-decostand(DivIndices[,11:16], method = "standardize", na.rm = T) %>%
  cbind(plot_field_id) %>%
  left_join(SiteList) %>%
  left_join(DivIndices[,1:10])

DivInd.stand2 <- DivInd.stand %>%
  gather(specMetric, specDiv.stand, CV.Casi:SV.Sasi) %>%
  separate(specMetric, c("specMetric",'sensor'))

# This function creates Spectral Diversity ~ Field Diversity plots
DivPlots <- function(df, metric, spectral){
  Plot <- ggplot(data = df %>% dplyr::filter(sensor == {{spectral}}),
                 aes(y = specDiv.stand, x = {{metric}})) +
    geom_point(aes(fill = specMetric, shape = specMetric)) +
    geom_smooth(aes(colour = specMetric, fill = specMetric),method = 'lm', fullrange = T, se = F) +
    scale_shape_manual(values = c(21,22,24),  name='Spectral Metric') +
    scale_colour_manual(values = c("black","#636363","#bdbdbd"), name = 'Spectral Metric') +
    scale_fill_manual(values = c("black","#636363","#bdbdbd"), name = 'Spectral Metric') +
    scale_y_continuous(expand=c(0,0), breaks=seq(-2, 4, by = 1), limits=c(-2.5,4)) +
    ylab("") +
    theme_classic() +
    theme(axis.title.x=element_markdown(face='bold', size=12),
          axis.text.x = element_text(size=10, colour = 'black'),
          title = element_text(face = 'bold', size=10),
          axis.text.y=element_text(size=10, colour = 'black'),
          axis.title.y =element_text(face='bold', size=12),
          legend.title=element_text(face='bold', size = 10),
          legend.position = 'none',
          #legend.box.background = element_rect(colour = "black"),
          legend.title.align=0.5)
  return(Plot)
  
}

DivIndPlots <- ggarrange(DivPlots(DivInd.stand2, q0, 'Casi') +
                    labs(title = 'Species Richness', x = '<sup>0</sup>*D*') +
                    scale_x_continuous(expand = c(0,0), limits = c(0,15), breaks = seq(0,15, by =5)),
                  DivPlots(DivInd.stand2, q1, 'Casi') +
                    labs(title = 'Exp. Shannon Index', x = '<sup>1</sup>*D*') +
                    scale_x_continuous(expand = c(0,0),limits = c(0,8), breaks = seq(0,8, by = 2)),
                  DivPlots(DivInd.stand2, q2, 'Casi') +
                    labs(title = 'Inv. Simpson Index', x = '<sup>2</sup>*D*') +
                    scale_x_continuous(expand = c(0,0), limits = c(0,6), breaks = seq(0,6,by = 2)),
                  DivPlots(DivInd.stand2, J, 'Casi') +
                    labs(title = 'Jaccards Evenness', x = "J\'") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = .25)),
                  DivPlots(DivInd.stand2, FRic, 'Casi') +
                    labs(title = 'Functional Richness', x = "F<sub>Ric</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,33), breaks= seq(0,30, by = 10)),
                  DivPlots(DivInd.stand2, FEve, 'Casi') +
                    labs(title = 'Functional Evenness', x = "F<sub>Eve</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by =0.25)),
                  DivPlots(DivInd.stand2, FDiv, 'Casi') +
                    labs(title = 'Functional Divergence', x = "F<sub>Div</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0.25,1), breaks = seq(0.25,1, by = 0.25)),
                  DivPlots(DivInd.stand2, FDis, 'Casi') +
                    labs(title = "Functional Dispersion", x ="F<sub>Dis</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                  DivPlots(DivInd.stand2, RaoQ, 'Casi') +
                    labs(title = 'Rao Entropy') +
                    scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = seq(0,12, by = 4)),
                  common.legend = TRUE, legend = 'bottom', labels = 'AUTO'
) 
DivIndPlots <- annotate_figure(DivIndPlots, left = text_grob("Standardized VNIR-Spectral Diversity", color = 'black', face = 'bold', size = 12, rot = 90))

DivIndPlots

ggsave(plot = DivIndPlots, filename = 'Outputs/Figures/DiversityCor_Casi.jpeg', width= 7.08661, height = 9.84252, dpi=600)

rm(DivInd.stand, DivInd.stand2, DivIndPlots)

# 2.0__Environmental Drivers of Community Properties ---------------------------------------------
#   Compare whether biodiversity dimensions are driven by the same environmental gradients - AT MONT MEGANTIC ONLY
#       - RDAs for composition and GAMs for diversity.

# But first, examine the correlation between environmental variables
EnvrCor <- cor(Envr[,8:13], use="complete.obs") #Do correlation, where we exclude Na values

corrplot(EnvrCor, method = 'number') # Roughness and slope are strongly correlated ... don't model roughness!

# Make sure environmental variables are dimensionally homogeneous (e.g., centered and standardized), so regression coefficients can be directly compared.
Envr.MMG <- Envr %>%
  dplyr::filter(site =='MtMeg-1')

Envr.MMG.stand <- decostand(Envr.MMG[,8:38], method = 'standardize') 

# 2.1__Composition  --------------------------------------------------------
# Redundancy analysis (RDA) is a direct gradient analysis method that summarises linear relationships between components of response variables that are explained by a set of explanatory variables.
# The total variance of the data set is partitioned into constrained and unconstrained variances,
#   -If constrained variance > unconstrained variance then analysis suggests that variation in response variables is accounted by explanatory variables.
#   -If unconstrained variance > constrained variance then the results should be interpreted with caution.
# Significance value for the overall RDA solution and individual RDA axes can be determined by permutation.

# This function fits RDAs and creates dataframes with different summary statistics related to:  
#   a) the overall model,
#   b) the canonical axes, 
#   c) the explanatory variables


# 2.1.1__Fit RDAs ---------------------------------------------------------

RDA <- function(df1, df2, site, dimension){
  rda <- rda(df1 ~ elevation + slope  + northness + eastness + TWI, data = df2)
  rda.summary <- summary(rda)
  rda.anova <- anova.cca(rda, permutations = 999)
  rda.anova.axes <- anova.cca(rda, by = 'axis', permutations = 999)
  rda.anova.mar <- anova(rda, by = 'mar', permuations = 999) 
  
  adj.R2 <- RsquareAdj(rda)$adj.r.squared
  fvalue <- rda.anova[1,3]
  pvalue <- rda.anova[1,4]
  constrained <- rda.summary[['constr.chi']] / rda.summary[['tot.chi']]
  unconstrained <- rda.summary[["unconst.chi"]] / rda.summary[['tot.chi']]
  
  rda.output <- data.frame(site = site,
                           dimension = dimension,
                           adj.R2 = adj.R2,
                           fvalue = fvalue,
                           pvalue = pvalue,
                           constrained = constrained,
                           unconstrained = unconstrained)
  
  rda.axes <- rda.anova.axes %>%
    rownames_to_column(var = 'axis') %>%
    mutate(site = {{site}}, .before = axis) %>%
    mutate(dimension = {{dimension}}, .before = axis)
  
  rda.explanvar <- rda.anova.mar %>% 
    rownames_to_column(var = 'explanvar') %>%
    mutate(site = {{site}}, .before = explanvar) %>%
    mutate(dimension = {{dimension}}, .before = explanvar)
  
  rda.correlations <- data.frame(rda[["CCA"]][['biplot']]) %>%
    rownames_to_column(var = 'explanvar') %>%
    mutate(site = {{site}}, .before = RDA1) %>%
    mutate(dimension = {{dimension}}, .before = RDA1)
  
  
  return(list(rda.output, rda.axes, rda.explanvar,rda.correlations, 
              rda))
}


# 2.1.1.1__Taxonomic Composition RDA ---------------------------------------------
TaxComp.sites <- TaxComp %>%
  left_join(SiteList) 

TaxComp.MMG <- TaxComp.sites %>%
  dplyr::filter(site == 'MtMeg-1')

TaxComp.MMG.hellinger <-sqrt(TaxComp.MMG[,3:32]) %>%
    dplyr::select(-c(10,14,19,25,28,30)) # This removes species that were not observed at Mont Megantic but were at Mont St Bruno

# Conduct RDAs 
Tax.MMG.rda <- RDA(TaxComp.MMG.hellinger, Envr.MMG.stand, 'MMG', 'Taxonomic')

# 2.1.1.2__Functional Composition RDA --------------------------------------------
FunComp.sites <- FunComp %>%
  left_join(SiteList)

FunComp.MMG <- FunComp.sites %>%
  dplyr::filter(site =='MtMeg-1')

FunComp.MMG.stand <- decostand(FunComp.MMG[,2:16], method = 'standardize') 

# Conduct RDAS
Fun.MMG.rda <- RDA(FunComp.MMG.stand, Envr.MMG.stand, 'MMG', 'Functional')

# Join output tables across dimensions
RDA.overall <- bind_rows(Tax.MMG.rda[[1]], Fun.MMG.rda[[1]])
RDA.axes <- bind_rows(Tax.MMG.rda[[2]], Fun.MMG.rda[[2]])
RDA.explanvar <- bind_rows(Tax.MMG.rda[[3]], Fun.MMG.rda[[3]])
RDA.correlations <- bind_rows(Tax.MMG.rda[[4]], Fun.MMG.rda[[4]])


# 2.1.1.3__VNIR-Spectral Composition RDA --------------------------------------------------
CasiComp.MMG <- CasiComp %>%
  dplyr::filter(Site == 'Megantic')

CasiComp.MMG.mean <- CasiComp.MMG %>% 
  select(c(6,20:248)) %>%  # Spectral reflectance is in columns 20:248
  group_by(plot_field_id) %>%
  summarise_all(mean)

Casi.MMG.rda <- RDA(CasiComp.MMG.mean[,-1], Envr.MMG.stand, 'MMG', 'Casi')

# Join output tables across dimensions
RDA.overall   <- bind_rows(RDA.overall, Casi.MMG.rda[[1]])
RDA.axes      <- bind_rows(RDA.axes, Casi.MMG.rda[[2]]) 
RDA.explanvar <- bind_rows(RDA.explanvar, Casi.MMG.rda[[3]])
RDA.correlations <- bind_rows(RDA.correlations, Casi.MMG.rda[[4]])

# 2.1.1.4__SWIR-Spectral Composition RDA --------------------------------------------------
SasiComp.MMG <- SasiComp %>%
  dplyr::filter(Site == 'Megantic')

SasiComp.MMG.mean <- SasiComp.MMG %>% 
  select(c(6,22:99)) %>%  # Spectral reflectance is in columns 20:248
  group_by(plot_field_id) %>%
  summarise_all(mean)

Sasi.MMG.rda <- RDA(SasiComp.MMG.mean[,-1], Envr.MMG.stand, 'MMG', 'Sasi')

# Join output tables across dimensions
RDA.overall <- bind_rows(RDA.overall, Sasi.MMG.rda[[1]])
RDA.axes      <- bind_rows(RDA.axes, Sasi.MMG.rda[[2]]) 
RDA.explanvar <- bind_rows(RDA.explanvar, Sasi.MMG.rda[[3]])
RDA.correlations <- bind_rows(RDA.correlations, Sasi.MMG.rda[[4]])

write.csv(RDA.overall, file = "Outputs/Statistics/RDA_Stats_Overall.csv")
write.csv(RDA.axes, file = "Outputs/Statistics/RDA_Stats_Axes.csv")
write.csv(RDA.explanvar, file = 'Outputs/Statistics/RDA_Stats_ExplanatoryVars.csv')
write.csv(RDA.correlations, file = 'Outputs/Statistics/RDA_Stats_Correlations.csv')


# 2.1.2__Plot RDAs ------------------------------------------------------
# Make RDA triplots for Taxonomic, Functional, and Casi for MMG
# Determine the percentage of conifer abundance in each plot - to aide the interpretation of trioplots
PercConifer.MMG <- PercConifer %>%
  filter(site == "MtMeg-1")

# This function fits the rda and returns statistics which will be used in the triplots
plotRDA1 <- function(df1, df2){
  rda <- rda(df1 ~ elevation + slope + northness + eastness + TWI, data = df2)
  rda.summary <- summary(rda)
  rda.anova <- anova(rda, permutations = 999)
  rda.anova.mar <- anova(rda, by = 'mar', permuations = 999) 
  
  adj.R2 <- RsquareAdj(rda)$adj.r.squared
  pvalue <- rda.anova[1,4]
  
  rda.plot <- ordiplot(rda, scaling = 2)
  # siteScores <- data.frame(rda.plot[['sites']]) %>% 
  #      cbind(PercCon = PercConifer.MMG$PercCon)
  spScores.rda <- data.frame(rda.plot[['species']])
  envrScores <- data.frame(rda.plot[['biplot']]) %>%
        rownames_to_column() %>%
        mutate(RDA1 = RDA1*attr(rda.plot[['biplot']],"arrow.mul"),
               RDA2 = RDA2*attr(rda.plot[['biplot']],"arrow.mul"))
    
  siteScores <- data.frame(rda.summary$constraints[,1:2]) %>%
    cbind(PercCon = PercConifer.MMG$PercCon)
  # spScores.rda <- data.frame(rda.summary$species[,1:2])
  # envrScores <- data.frame(rda.summary$biplot[,1:2]) 
  
  EnvFit <- envfit(rda, env = df1) 
  spScores <- data.frame(r = EnvFit[["vectors"]][["r"]],
                         pval = EnvFit[["vectors"]][["pvals"]]) %>%
                cbind(spScores.rda) %>%
                arrange(desc(r)) %>%
                slice(1:5)
  
  RDA1 <- rda.summary$cont[[1]][2,1] / rda.summary$cont$importance[3,5] * adj.R2
  RDA2 <- rda.summary$cont$importance[2,2] / rda.summary$cont$importance[3,5] * adj.R2
  
  return(list(adj.R2, pvalue, siteScores, spScores, envrScores, RDA1, RDA2) %>%
           set_names(c('adj.R2', 'pvalue', 'siteScores', 'spScores', 'envrScores', 'RDA1', 'RDA2')))
}

# This function plots the triplots
plotRDA2 <- function(plotRDA1.output){
  biplot <- ggplot()+
                geom_point(data = plotRDA1.output[['siteScores']],
                           aes(x = RDA1, y = RDA2, fill = PercCon), 
                           size = 2, pch = 21) +
               # scale_fill_brewer(palette = 'RdYlGn') +
                scale_fill_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
    
                geom_segment(data = plotRDA1.output[['spScores']],
                              aes(x=0, y=0, xend=RDA1, yend=RDA2),
                              arrow = arrow(type = 'closed', length = unit(0.1,"cm")),
                              linetype = 'dashed', 
                              colour = 'grey30') +
                # geom_text(data = plotRDA1.output[['spScores']],
                #               aes(x  = RDA1, y = RDA2, label = row.names(plotRDA1.output[['spScores']])),
                #               size=4) +

                geom_segment(data = plotRDA1.output[['envrScores']],
                             aes(x=0, y=0, xend=RDA1, yend=RDA2),
                             arrow = arrow(type = 'closed', length = unit(0.1,"cm")), 
                             colour = 'black') +
                # geom_text_repel(data = plotRDA1.output[['envrScores']],
                #           aes(x  = RDA1, y = RDA2, label = row.names(plotRDA1.output[['envrScores']])), 
                #           size=4) +
                # annotate('text', x = plotRDA1.output[['envrScores']][1,1], y =plotRDA1.output[['envrScores']][1,2], label = 'Elevation') +
                # annotate('text', x = plotRDA1.output[['envrScores']][2,1], y =plotRDA1.output[['envrScores']][2,2], label = 'Slope') +
                # annotate('text', x = plotRDA1.output[['envrScores']][3,1], y =plotRDA1.output[['envrScores']][3,2], label = 'Northness') +
                # annotate('text', x = plotRDA1.output[['envrScores']][4,1], y =plotRDA1.output[['envrScores']][4,2], label = 'Eastness') +
                # annotate('text', x = plotRDA1.output[['envrScores']][5,1], y =plotRDA1.output[['envrScores']][5,2], label = 'TWI') +
                # 
                geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
                geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
                labs(x=paste("RDA 1 (", format(100 *plotRDA1.output[['RDA1']], digits=4), "%)", sep=""),
                     y=paste("RDA 2 (", format(100 *plotRDA1.output[['RDA2']], digits=4), "%)", sep="")) +
                #annotate('text', x=0.4, y=-0.5, label= paste("adj.italic(R) ^2 == (", format(plotRDA1.output[['adj.R2']], digits =3), ")")) + 
                theme_bw() +
                theme(panel.grid = element_blank(), 
                      axis.text.x = element_text(size=8, color='black'),
                      axis.text.y = element_text(size=8, color='black'),
                      axis.title.y = element_text(face='bold', size=10),
                      axis.title.x = element_text(face='bold', size=10),
                      legend.title = element_text(face='bold', size = 10),
                      legend.text = element_text(face = 'bold', size= 10, color ='black'),
                      legend.position = 'none') 
  return(biplot)              
}

# Create triplot, use annotate to position environ "species" arrows nicely - NOTE: I ended up using an illustrated to add the annotations
windows(width = 7.08661, height = 9.44882)
RDAbiplots <- ggarrange(plotRDA2(plotRDA1(TaxComp.MMG.hellinger, Envr.MMG.stand)) +
                          labs(title = 'Taxonomic') + 
                          theme(plot.title = element_text(face = "bold")) +
                          # annotate('text', x=-0.8, y=-1, label="adj.~italic(R) ^2 == 0.27", parse = TRUE, size = 4) +
                          guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
                          scale_x_reverse(limits = c(0.8, -0.8)) + 
                          scale_y_continuous(limits = c(-0.8, 0.8)) ,
                          # annotate('text', x = (-0.7), y = (0.7), label="Elevation", size = 4, fontface = 2) +
                          # annotate('text', x = -0.18562583, y = 0.77, label="Slope", size = 4, fontface = 2) +
                          # annotate('text', x = -0.55, y = -0.7, label="Northness", size = 4, fontface = 2) +
                          # annotate('text', x = -0.17, y = -0.08, label="Eastness", size = 4) +
                          # annotate('text', x = 0.15975848, y = -0.68, label="TWI", size = 4) +
                          # annotate('text', x = 0.85, y = 0.26314013, label="ACSA", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.72024377, y = 0.26, label="ABBA", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.48, y = -0.075, label="BEPA", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.05, y = -0.28, label="PIRU", size = 4, colour = 'grey30') +
                          # annotate('text', x = 0.19, y = -0.25, label="BEAL", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.11388143, y = 0.01079594, label = "SODE", size = 4, colour = 'grey30') +
                          # annotate('text', x = 0.12736090, y = 0.01649647, label = "FAGR", size = 4, colour = 'grey30')  +
                          # annotate('text', x = 0.21645820, y = -0.21741349, label = "ACRU", size = 4, colour = 'grey30') +
                          # annotate('text', x = 0.05434935, y = -0.10870284, label = "POTR", size = 4, colour = 'grey30') ,
                        plotRDA2(plotRDA1(FunComp.MMG.stand, Envr.MMG.stand)) + 
                          labs(title = 'Functional') + 
                          theme(plot.title = element_text(face = "bold")) +
                          scale_y_continuous(limits = c(-2.2, 2)) ,
                          # scale_x_continuous(limits = c(-1.5, 2)) ,
                          # annotate('text', x=1.8, y=-2.8, label="adj.~italic(R) ^2 == 0.28", parse = TRUE, size = 4) +
                          # annotate('text', x = 0.78523944, y = 0.68, label="Elevation", size = 4, fontface = 2) +
                          # annotate('text', x = 0.16866749, y = 0.72, label="Slope", size = 4, fontface = 2) +
                          # annotate('text', x = 0.3, y = -0.8, label="Northness", size = 4, fontface = 2) +
                          # annotate('text', x = 0.35, y = -0.1762961, label="Eastness", size = 4) +
                          # annotate('text', x = -0.25, y = -0.88, label="TWI", size = 4) +
                          # annotate('text', x = -1.05, y = 0.1, label="SLA", size = 4, colour = 'grey30') +
                          # annotate('text', x = 1.05, y = -0.09026286, label="EWT", size = 4, colour = 'grey30') +
                          # annotate('text', x = -1.1, y = -0.09, label="Carot", size = 4, colour = 'grey30') +
                          # annotate('text', x = -1, y = -0.36, label="Chl b", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.93, y = -0.52, label="Chl a", size = 4, colour = 'grey30'), 
                        plotRDA2(plotRDA1(CasiComp.MMG.mean[,-1], Envr.MMG.stand)) + 
                          labs(title ='VNIR-Spectral') + 
                          theme(plot.title = element_text(face = "bold")) +
                          # annotate('text', x=1.35, y=2, label="adj.~italic(R) ^2 == 0.44", parse = TRUE, size = 4) +
                          scale_y_reverse(limits = c(1.2, -1.2)) + 
                          scale_x_continuous(limits = c(-0.8, 1.7), breaks = c(-0.5, 0, 0.5, 1, 1.5)) , #+
                          # annotate('text', x = 0.86, y = -0.34, label="Elevation", size = 4, fontface = 2) +
                          # annotate('text', x = 0.36704932, y = -0.76, label="Slope", size = 4, fontface = 2) +
                          # annotate('text', x = -0.16381222, y = 0.68, label="Northness", size = 4) +
                          # annotate('text', x = 0.18, y = -0.15, label="Eastness", size = 4) +
                          # annotate('text', x = -0.5, y = 0.75, label="TWI", size = 4) +
                          # annotate('text', x = -0.44, y = 0.15, label="466 nm", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.41, y = 0.05, label="468 nm", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.23, y = -0.16, label="523 nm", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.25, y = -0.05, label="521 nm", size = 4, colour = 'grey30') +
                          # annotate('text', x = -0.23, y = -0.27, label="573 nm", size = 4, colour = 'grey30'), 
                        nrow = 3, ncol = 1,
                        common.legend = TRUE, legend = 'bottom',
                        labels =  'AUTO'
                        ) 
RDAtriplots

ggsave(plot = RDAtriplots, filename = 'Outputs/Figures/RDA_triplots.jpeg', width=3.346, height=9.44882, dpi=600)

rm(CasiComp.MMG, CasiComp.MMG.mean, Casi.MMG.rda, EnvrCor, Fun.MMG.rda, FunComp.MMG, FunComp.MMG.stand, FunComp.sites, RDA.axes, 
   RDA.correlations, RDA.explanvar, RDA.overall, RDAbiplots, SasiComp.MMG, SasiComp.MMG.mean, Sasi.MMG.rda, TaxComp.MMG.hellinger, TaxComp.sites, Tax.MMG.rda)


# 2.2__Diversity -----------------------------------------------------------
# Generalized Additive Models (GAMs) is a GLM in which the linear response variable depends linearly on unknown smooth functions of some predictor variables. 
# The relationship between individual predictors and the response variable follow smooth patterns that can be linear or non-linear. There are two pros:
#     i)  Flexible predictor functions can uncover hidden patterns in data
#     ii) Regularization of predictor functions helps avoid overfitting. 
# There are two ways of estimating the smoothing parameter:
#     i)  Generalized cross-validation (GCV)  
#     ii) Restricted maximum likelihood (REML) - Here, I will use REML because GCV is prone to under-smoothing
# Note: package mgcv uses penalized regression splines and by default defines basis functions 


# 2.2.1__Fit GAMs ---------------------------------------------------------

# Subset data for just Mont Megantic and merge the diversity dataframe with the standardized environmental variable dataframe
DivIndices.MMG <- DivIndices %>%
  left_join(SiteList) %>%
  dplyr::filter(site == 'MtMeg-1') %>%
  cbind(Envr.MMG.stand) %>%
  left_join(PercConifer.MMG)


# 2.2.1.1__Taxonomic Diversity GAM ----------------------------------------
# Taxonomic Diversity(q1)
GAM.q1  <- gam(q1 ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.q1)

Check.GAM.q1 <- getViz(GAM.q1) 
check(Check.GAM.q1,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.q1)

concurvity(GAM.q1, full = F)

plot(GAM.q1, residuals = T, pch = 19, cex=.3)


# 2.2.1.2__Functional Diversity GAM ---------------------------------------
# Functional Diversity (FDis)
GAM.FDis  <- gam(FDis ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.FDis)

Check.GAM.FDis <- getViz(GAM.FDis) 

check(Check.GAM.FDis,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.FDis)

plot.gam(GAM.FDis, residuals = T, pch = 19, cex=.3)


# 2.2.1.3_VNIR-Spectral Diversity GAM -------------------------------------
# CASI Diversity (SV.Casi)
GAM.SV.Casi  <- gam(SV.Casi ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.SV.Casi)

Check.GAM.SV.Casi <- getViz(GAM.SV.Casi) 

check(Check.GAM.SV.Casi,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.SV.Casi)

plot(GAM.SV.Casi, residuals = T, pch = 19, cex=.3)


# 2.2.1.4__SWIR-Spectral Diversity GAM ------------------------------------
# SASI Diversity (SV.Sasi)
GAM.SV.Sasi  <- gam(SV.Sasi ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.SV.Sasi)

Check.GAM.SV.Sasi <- mgcViz::getViz(GAM.SV.Sasi) 

check(Check.GAM.SV.Sasi,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.SV.Sasi)

plot(GAM.SV.Sasi, residuals = T, pch = 19, cex=.3)

# GAMextract - a function for quickly extracting GAM model outputs
GAMextract <- function(df, metric, site, dimension){
  df2 <- df %>%
    dplyr::select(c(metric, 18:23)) %>%
    #select(c(metric, 17:22)) %>%
    rename(index = metric)
  
  full <- gam(index ~ s(elevation) + TWI + northness + eastness + slope, method = 'REML', data = df2)
  elevation <- gam(index ~ TWI + northness + eastness + slope, method = 'REML', data = df2)
  TWI <- gam(index ~ s(elevation) + northness + eastness + slope, sp = c(full$sp[1]), method = 'REML', data = df2)
  northness <- gam(index ~ s(elevation) + TWI + eastness + slope, sp = c(full$sp[1]), method = 'REML', data = df2)
  eastness <- gam(index ~ s(elevation) + TWI + northness + slope, sp = c(full$sp[1]), method = 'REML', data = df2)
  slope <- gam(index ~ s(elevation) + TWI + northness + eastness, sp = c(full$sp[1]), method = 'REML', data = df2)
  null <- gam(index ~ 1, data = df2)

  devElev <- (elevation$deviance - full$deviance) / null$deviance
  devTWI  <- (TWI$deviance - full$deviance) / null$deviance
  devNorthness <- (northness$deviance - full$deviance) / null$deviance
  devEastness <- (eastness$deviance - full$deviance) / null$deviance
  devSlope <- (slope$deviance - full$deviance) / null$deviance

  modelsummary <- summary(full)
  adj.R2 <- modelsummary$r.sq
  dev.expl <- modelsummary$dev.expl
  AIC <- AIC(full)
  GCV <- full[["gcv.ubre"]][["REML"]]
  
  GAMoutput <- data.frame(site = site,
                           dimension = dimension,
                           index = metric,
                           adj.R2 = adj.R2,
                           dev.expl = dev.expl,
                           AIC = AIC,
                           GCV = GCV)
  
  GAMexplanvar <- data.frame(metric = metric,
                             explanvar = c('Elevation', 'TWI', 'Northness', 'Eastness', 'Slope'),
                             deviance = c(devElev, devTWI, devNorthness, devEastness, devSlope)) %>%
    mutate(proportion = deviance/sum(deviance) * 100)

  return(list(GAMoutput, GAMexplanvar))
}

# Fit GAMs for q1, FDis, and SV
Tax.MMG.gam <- GAMextract(DivIndices.MMG, 'q1', 'MMG', 'Taxonomic')
Fun.MMG.gam <- GAMextract(DivIndices.MMG, 'FDis', 'MMG', 'Functional')
Casi.MMG.gam <- GAMextract(DivIndices.MMG, 'SV.Casi', 'MMG', 'Casi')
Sasi.MMG.gam <- GAMextract(DivIndices.MMG, 'SV.Sasi', 'MMG', 'Sasi')

GAM_overall <- bind_rows(Tax.MMG.gam[[1]], Fun.MMG.gam[[1]], Casi.MMG.gam[[1]], Sasi.MMG.gam[[1]])
GAM_explanvar <- bind_rows(Tax.MMG.gam[[2]], Fun.MMG.gam[[2]], Casi.MMG.gam[[2]], Sasi.MMG.gam[[2]])

write.csv(GAM_overall, file = "Outputs/Statistics/GAM_Stats_Overall.csv")
write.csv(GAM_explanvar, file = "Outputs/Statistics/GAM_Stats_Explanvar.csv")

# 2.2.2__Plot GAMs ---------------------------------------------------------

# Create a dataframe with unstandardized environmental variables
DivIndices.MMG.raw <- DivIndices %>%
  left_join(SiteList) %>%
  dplyr::filter(site == 'MtMeg-1') %>%
  merge(Envr.MMG, by.y = 'field_plot_id', by.x = 'plot_field_id') %>%
  left_join(PercConifer.MMG)

# This function plots the relationship with elevation (ie. the smooth variable in the GAMs)
plotGAM <- function(metric, envrvar, test, df, df3){
  df2 <- df %>%
    dplyr::select(c({{metric}},18:23)) %>%
    rename(index = {{metric}})
  
  full <- gam(index ~ s(elevation) + TWI + northness + eastness + slope, method = 'REML', data = df2)
  
  smooth <- ifelse(test == 'elevation', full$sp[1],
                   ifelse(test == 'TWI', full$sp[2],
                          ifelse(test == 'northness', full$sp[3],
                                 ifelse(test == 'eastness', full$sp[4], full$sp[5]))))
  
  Plot<-ggplot(df3, aes(x={{envrvar}}, y={{metric}}))+
    stat_smooth(method=gam, formula=y~s(x, sp = smooth), se=T,colour='black') +
    #stat_smooth(method = lm, se=T, colour='black') +
    geom_point(aes(fill = PercCon), size = 2, pch = 21) +
    scale_fill_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
    #scale_y_continuous(expand=c(0,0)) + 
    #scale_x_continuous(expand=c(0,0)) + 
    theme_classic() +
    theme(axis.title.x=element_text(face='bold', size=12, colour = 'black'),
          axis.text.x = element_text(size=10, colour = 'black'),
          axis.text.y=element_text(size=10, colour = 'black'),
          axis.title.y =element_markdown(face='bold', size=12, colour = 'black'),
          legend.text = element_text(face = 'bold', size = 10, colour = 'black'), 
          legend.title=element_text(face='bold', size = 10),
          plot.title = element_text(face = 'bold'),
          legend.title.align=0.5,
          legend.justification = 'center'
          )
}

# This function plots the relationship with the linear fitted variables
plotGAM2 <- function(metric, envrvar, test, df, df3){
  df2 <- df %>%
    dplyr::select(c({{metric}},18:23)) %>%
    rename(index = {{metric}})
  
  full <- gam(index ~ s(elevation) + TWI + northness + eastness + slope, method = 'REML', data = df2)
  
  smooth <- ifelse(test == 'elevation', full$sp[1],
                   ifelse(test == 'TWI', full$sp[2],
                          ifelse(test == 'northness', full$sp[3],
                                 ifelse(test == 'eastness', full$sp[4], full$sp[5]))))

  Plot<-ggplot(df3, aes(x={{envrvar}}, y={{metric}}))+
    #stat_smooth(method=gam, formula=y~s(x), se=T,colour='black') +
    stat_smooth(method = lm, se=T, colour='black') +
    geom_point(aes(fill = PercCon), size = 3, pch = 21) +
    scale_fill_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
    #scale_y_continuous(expand=c(0,0)) + 
    #scale_x_continuous(expand=c(0,0)) + 
    theme_classic() +
    theme(axis.title.x=element_text(face='bold', size=12, colour = 'black'),
          axis.text.x = element_text(size=10, colour = 'black'),
          axis.text.y=element_text(size=10, colour = 'black'),
          axis.title.y =element_markdown(face='bold', size=12, colour = 'black'),
          legend.title=element_text(face='bold', size = 10),
          plot.title = element_text(face = 'bold'),
          legend.title.align=0.5,
          legend.justification = 'center'
    )
}

GAMPlot <- ggarrange(plotGAM(q1, elevation, 'elevation', DivIndices.MMG, DivIndices.MMG.raw) +
                       labs(title = 'A. Taxonomic', x='', y='Exp. Shannon Index (<sup>1</sup>*D*)') +
                       scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)) ,
                     plotGAM2(q1, northness, 'northness', DivIndices.MMG, DivIndices.MMG.raw) +
                       labs(title = '', x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)) ,
                     plotGAM2(q1, eastness, 'eastness', DivIndices.MMG, DivIndices.MMG.raw) +
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)),
                     plotGAM2(q1, slope, 'slope', DivIndices.MMG, DivIndices.MMG.raw) +
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(0,30), breaks = seq(0,30, by = 5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)),
                     plotGAM2(q1, TWI, 'TWI', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '', x='', y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(2,12), breaks = seq(0,12, by = 2)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)),

                     plotGAM(FDis, elevation, 'elevation', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = 'B. Functional', x='', y='Functional Dispersion (F<sub>Dis</sub>)') +
                        scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                        scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM2(FDis, northness, 'northness', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM2(FDis, eastness, 'eastness', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM2(FDis, slope, 'slope', DivIndices.MMG, DivIndices.MMG.raw) +
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(0,30), breaks = seq(0,30, by = 5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM2(FDis, TWI, 'TWI', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(2,12), breaks = seq(0,12, by = 2)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),

                     plotGAM(SV.Casi, elevation, 'elevation', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = 'C. VNIR-Spectral', x='Elevation\n(m a.s.l.)', y='Spectral Variance (SV)') +
                       scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM2(SV.Casi, northness, 'northness', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '',x='Northness\n(cos(Aspect))',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM2(SV.Casi, eastness, 'eastness', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '',x='Eastness\n(sin(Aspect))',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM2(SV.Casi, slope, 'slope', DivIndices.MMG, DivIndices.MMG.raw) +
                       labs(title = '',x='Slope\n(\u00B0)',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(0,30), breaks = seq(0,30, by = 5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM2(SV.Casi, TWI, 'TWI', DivIndices.MMG, DivIndices.MMG.raw)+
                       labs(title = '',x='Topographic Wetness Index\n',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(2,12), breaks = seq(0,12, by = 2)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),

                     nrow = 3, ncol = 5,
                     common.legend = T, legend = 'bottom'
                     )

GAMPlot

ggsave(plot = GAMPlot, filename = 'Outputs/Figures/GAM_plots.jpeg', width=14, height=11, dpi=300)


GAMElevation <- ggarrange(plotGAM(q1, elevation, 'elevation', DivIndices.MMG, DivIndices.MMG.raw) +
                             labs(title = 'Taxonomic', x='Elevation (m a.s.l.)', y='Exp. Shannon Index (<sup>1</sup>*D*)') +
                             scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                             scale_y_continuous(expand = c(0,0), limits = c(0,7), breaks = seq(0,7, by = 1)), 
                          
                          plotGAM(FDis, elevation, 'elevation', DivIndices.MMG, DivIndices.MMG.raw)+
                            labs(title = 'Functional',  x='Elevation (m a.s.l.)', y='Functional Dispersion (F<sub>Dis</sub>)') +
                            scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                            scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                          
                          plotGAM(SV.Casi, elevation, 'elevation', DivIndices.MMG, DivIndices.MMG.raw)+
                            labs(title = 'VNIR-Spectral', x='Elevation (m a.s.l.)', y='Spectral Variance (SV)') +
                            scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) ,
                            #scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                          
                          nrow = 3, ncol = 1,
                          common.legend = T, legend = 'bottom',
                          labels =  'AUTO')

GAMElevation

ggsave(plot = GAMElevation, filename = 'Outputs/Figures/GAM_plots_elevation2.jpeg', width=3.346, height=9.44882, dpi=600)

rm(Casi.MMG.gam, Check.GAM.q1, Check.GAM.FDis, Check.GAM.SV.Casi, Check.GAM.SV.Sasi, Fun.MMG.gam, GAM_explanvar, GAM_overall, GAM.FDis,
   GAM.q1, GAM.SV.Casi, GAM.SV.Sasi, GAMElevation, GAMPlot, Sasi.MMG.gam, Tax.MMG.gam)

# 3.0__Extra Figures ------------------------------------------------------

# 3.1__Study Site Figure -------------------------------------------------------

# Plot the Spectral Composition for Plot 1 
CasiPlot1 <- CasiComp %>%
  dplyr::filter(plot_field_id == 'Plot01') %>%
  dplyr::select(-c(1:19)) %>%
  gather(Wavelength, BandDepth) %>%
  mutate(Wavelength = gsub("BD_Band", "", Wavelength),
         Wavelength = gsub("nm", "", Wavelength)) %>%
  group_by(Wavelength) %>%
  summarise(BDavg = mean(BandDepth),
            #BDsd = sd(BandDepth),
            BDmax = max(BandDepth),
            BDmin = min(BandDepth)) %>%
  mutate(Sensor = 'Casi',
         Wavelength = as.numeric(Wavelength))

SasiPlot1 <- SasiComp %>%
  dplyr::filter(plot_field_id == 'Plot01') %>%
  dplyr::select(-c(1:21)) %>%
  gather(Wavelength, BandDepth) %>%
  mutate(Wavelength = gsub("BD_Band", "", Wavelength),
         Wavelength = gsub("nm", "", Wavelength)) %>%
  dplyr::filter(Wavelength > 1059) %>%
  group_by(Wavelength) %>%
  summarise(BDavg = mean(BandDepth),
            #BDsd = sd(BandDepth),
            BDmax = max(BandDepth),
            BDmin = min(BandDepth)) %>%
  mutate(Sensor = 'Sasi',
          Wavelength = as.numeric(Wavelength))

NullWavelengths<- data.frame(Wavelength = c(899:957,1345:1460, 1790:1950)) %>%
                      mutate(Sensor = ifelse(Wavelength %in% (899:957), 'Casi', 'Sasi'))

SpectraPlot1 <- full_join(NullWavelengths, rbind(CasiPlot1, SasiPlot1))

SpectrumPlot <- ggplot(SpectraPlot1) +
                  #annotate("rect", xmin = 430, xmax = 1059, ymin = -0.08, ymax = 1, alpha = .3, color="grey40", fill= 'grey40', linetype = 'dashed') +
                  #annotate("rect", xmin = 972, xmax = 2435, ymin = -0.08, ymax = 1, alpha = .3, color="grey40", fill= 'grey20', linetype = 'dashed') +
                  #geom_mark_rect(aes(x=Wavelength, y=BDavg, fill = Sensor)) +
                  #scale_fill_manual(values = c('#9ccb66','#e9e29c')) +
                  #annotate("segment", x = 454, xend = 1059, y = -0.07, yend = -0.07, cex = 0.5) +
                  #annotate("text", x = 756.5, y = -0.06, vjust = 0, label = 'CASI', fontface =2, size = 5) +
                  #annotate("segment", x = 972, xend = 2412, y = -0.08, yend = -0.08, cex = 0.5) +
                  #annotate("text", x = 1692, y = -0.06, vjust = 0, label = 'SASI', fontface =2, size =5) +
                  geom_line(aes(x=Wavelength, y=BDavg, colour=Sensor), cex=1.1) +
                  geom_ribbon(aes(x= Wavelength, ymin=BDmin, ymax=BDmax, fill=Sensor), alpha=0.2) +
                  scale_colour_manual(values = c('#506969','#6e5353')) +
                  scale_fill_manual(values = c('#749898','rosybrown4')) +
                  scale_y_continuous(limits = c(-0.1,1), expand = c(0,0)) +
                  scale_x_continuous(limits = c(300,2500), expand = c(0,0), breaks = seq(400,2400,400)) +
                  theme_classic() +
                  xlab("Wavelength (nm)") +
                  ylab("Normalized Reflectance") +
                  theme(axis.title.x=element_text(face='bold', size=14),
                        axis.text.x = element_text(face='bold', size=12),
                        axis.text.y=element_text(face='bold', size=12),
                        axis.title.y =element_text(face='bold', size=14),
                        legend.title=element_text(face='bold', size = 12))
SpectrumPlot

ggsave(plot = SpectrumPlot, filename = 'Outputs/Figures/Spectra_Plot01.jpeg', width=5, height=4, dpi=300)

# Plot species composition for Plot 1
# library(sf)

# PointClouds<- ggplot() +
#                 geom_sf(data = TreeCircles %>% dplyr::filter(plot_field_id == 'Plot06_Carya'), alpha=0.3, color = 'black') +
#                 facet_wrap(~ Sensor) +
#                 geom_sf(data = CASIPointCloud %>% dplyr::filter(plot_field_id == 'Plot06_Carya') %>% mutate(Sensor = 'CASI'), aes(fill = Shadow), pch=21) +
#                 geom_sf(data = SASIPointCloud %>% dplyr::filter(plot_field_id == 'Plot06_Carya') %>% mutate(Sensor = 'SASI', Shadow = ShadowCASI2), aes(fill = Shadow), color = 'black', pch=21) +
#                 scale_fill_manual(values = c('grey80', 'black'), labels = c('Un-Masked', 'Masked')) +
#                 theme_bw() +
#                 theme(legend.position = 'none',
#                       strip.text.x = element_text(face='bold', size=12),
#                       axis.title.x=element_text(face='bold', size=12),
#                       axis.text.x = element_text(face='bold', size=10),
#                       axis.text.y=element_text(face='bold', size=10),
#                       axis.title.y =element_text(face='bold', size=12))
# 
# ggsave(plot = PointClouds, filename = 'Outputs/Figures/PointCloud_Plot06_Carya.png')



# 3.2__Northness, Elevation, and Red Spruce -------------------------------
test <- full_join(PercConifer.MMG, Envr.MMG) 

Northness <- full_join(TaxComp.MMG, PercConifer.MMG) %>%
  merge(Envr.MMG, by.x = 'plot_field_id', by.y = 'field_plot_id') %>%
  mutate(elevation = round(elevation, 0))

Northness_plot <- ggplot(Northness, aes(x = northness, y = elevation, fill = PercCon, label = elevation))+
  geom_point(size = 2, pch = 21) +
  #geom_text(hjust = 0.5 , vjust = -0.5) +
  scale_fill_gradient2(low = 'white', mid = 'black', high = 'white', midpoint = 0.5, name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
  theme_classic() +
  ylab("Elevation (m a.s.l.) ")+
  xlab("Northness\n (cos(aspect))")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
  labs(title = "") +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.title =element_text(face='bold', size=10),
        legend.text = element_text(face='bold', size=10),
        legend.title = element_text(face='bold', size=10),
        legend.position = 'bottom',
        legend.title.align=0.5) 
Northness_plot

Northness_plot2 <- ggplot(Northness, aes(x = northness, y = Picea.rubens.Sargent, label = elevation, fill = elevation)) +
  geom_point(size = 2, pch = 21) +
  #geom_text(hjust = 0.5 , vjust = -0.5) +
  scale_fill_gradient2(low = 'white', mid = 'black', high = 'white', midpoint = median(Envr.MMG$elevation), name = "Elevation (m a.s.l)") +
  theme_classic() +
  ylab("Relative Abundance\n of Picea rubens (%) ")+
  xlab("Northness\n (cos(aspect))")+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
  labs(title = "") +
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=10, colour = 'black'),
        axis.text.y=element_text(size=10, colour = 'black'),
        axis.title =element_text(face='bold', size=10),
        legend.text = element_text(face='bold', size=10),
        legend.title = element_text(face='bold', size=10),
        legend.position = 'bottom',
        legend.title.align=0.5) 
Northness_plot2

windows(width = 3.334646, height = 7)
Northness_plot3 <- ggarrange(Northness_plot,
                        Northness_plot2,
                        nrow = 2, ncol = 1,
                        common.legend = F, labels = "AUTO")
Northness_plot3

ggsave(plot = Northness_plot3, filename ='Outputs/Figures/Northness.jpeg', width=3.334646, height=7, dpi=600)


# 3.3__Wavelength Contributions to Spec Div -------------------------
#   Visualize wavelengths contributing to SV

# This function calculates the variance in normalized reflectance at each band per plot
calc.wavVar <- function(df, grouping_var, .x){
  
  points <- df %>%
    group_by({{grouping_var}}) %>%
    summarise(points = n())
  
  SV.df  <- df %>% group_by({{grouping_var}}) %>%
    summarise_all(~sum((.x - mean(.x))^2)) %>%
    rowwise({{grouping_var}}) %>%
    left_join(points) %>%
    summarise_all(~sum(.x / (points - 1)))
  
  return(SV.df)
}

Casi.SV2 <- CasiComp %>%
  dplyr::select(c(6,20:248)) %>%
  calc.wavVar(plot_field_id, c(2:230))

CasiSVPlot <- Casi.SV2 %>%
  gather(Wavelength, Variance, c(2:230)) %>%
  mutate(Wavelength = gsub("BD_Band", "", Wavelength),
         Wavelength = gsub("nm", "", Wavelength)) %>%
  group_by(Wavelength) %>%
  mutate(Varavg = mean(Variance)) %>%
  mutate(Wavelength = as.numeric(Wavelength)) %>%
  left_join(PercConifer.MMG) #%>%
#mutate(PercCon = as.character(PercCon))

Sasi.SV2 <- SasiComp %>%
  dplyr::select(c(6,22:99)) %>%
  calc.wavVar(plot_field_id, c(2:79))

SasiSVPlot <- Sasi.SV2 %>%
  gather(Wavelength, Variance, c(2:79)) %>%
  mutate(Wavelength = gsub("BD_Band", "", Wavelength),
         Wavelength = gsub("nm", "", Wavelength)) %>%
  group_by(Wavelength) %>%
  mutate(Varavg = mean(Variance)) %>%
  mutate(Wavelength = as.numeric(Wavelength)) %>%
  left_join(PercConifer.MMG)

windows(width = 3.35, height = 2)
CasiWav <- ggplot(CasiSVPlot) +
  geom_line(aes(x=Wavelength, y=Variance, group = plot_field_id, colour = PercCon)) +
  scale_colour_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
  #scale_colour_gradient2(low = 'white', mid = 'black', high = 'white', midpoint = 0.5, name = "Percent\n Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100 %')) +
  geom_line(aes(x=Wavelength, y=Varavg), colour = 'black', linewidth = 1.1) +
  annotate("rect", xmin = 899, xmax = 957, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
  ggbreak::scale_y_cut(breaks = 0.021, scales = 0.5, space  = .1) +
  scale_y_continuous(limits = c(-0.001,0.15), expand = c(0,0), breaks =c(0.0, 0.005, 0.01, 0.015, 0.02, .05, .10,.15)) +
  scale_x_continuous(limits = c(454,1059), breaks = seq(500,1000,100)) +
  theme_classic() +
  xlab("Wavelength (nm)") +
  ylab("Variance") +
  #labs(tag = "A") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(face='bold', size=12),
        axis.title.x = element_text(face='bold', size=12),
        legend.title = element_text(face='bold', size = 10),
        plot.tag = element_text(face = 'bold'),
        legend.position = 'none') 
CasiWav
ggsave(plot = CasiWav, filename = 'Outputs/Figures/Wavelenght_variance_VNIR.jpeg', width=3.334646, height=2, dpi=600)

SasiWav <- ggplot(SasiSVPlot) +
  geom_line(aes(x=Wavelength, y=Variance, group = plot_field_id, colour = PercCon)) +
  scale_colour_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
  #scale_colour_gradient2(low = 'white', mid = 'black', high = 'white', midpoint = 0.5, name = "Percent\n Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100 %')) +
  geom_line(aes(x=Wavelength, y=Varavg), colour = 'black', linewidth = 1.1) +
  #annotate("rect", xmin = 899, xmax = 957, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
  annotate("rect", xmin = 1345, xmax = 1460, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
  annotate("rect", xmin = 1790, xmax = 1950, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
  ggbreak::scale_y_cut(breaks = 0.021, scales = 0.5, space  = .1) +
  scale_y_continuous(limits = c(-0.001,0.15), expand = c(0,0), breaks =c(0.0, 0.005, 0.01, 0.015, 0.02, .05, .10,.15)) +
  scale_x_continuous(limits = c(1002,2412), breaks = seq(1000,2400,400)) +
  theme_classic() +
  xlab("Wavelength (nm)") +
  ylab("Variance") +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size=10, color='black'),
        axis.title.y = element_text(face='bold', size=12),
        axis.title.x = element_text(face='bold', size=12),
        legend.title = element_text(face='bold', size = 10),
        plot.tag = element_text(face = 'bold'),
        legend.position = 'none') 
SasiWav
ggsave(plot = SasiWav, filename = 'Outputs/Figures/Wavelenght_variance_SWIR.jpeg', width=3.334646, height=2, dpi=600)


# WavVar <- ggarrange(ggplot(CasiSDPlot) +
#                       labs(title = 'VNIR-Spectral') +
#                       guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
#                       geom_line(aes(x=Wavelength, y=Variance, group = plot_field_id, colour = PercCon)) +
#                       scale_colour_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
#                       geom_line(aes(x=Wavelength, y=Varavg), colour = 'black', linewidth = 1.1) +
#                       annotate("rect", xmin = 899, xmax = 957, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
#                       #ggbreak::scale_y_cut(breaks = 0.021, scales = 0.5, space  = .1) +
#                       scale_y_continuous(limits = c(-0.001,0.15), expand = c(0,0), breaks = seq(0.0, 0.15, 0.05)) +
#                       scale_x_continuous(limits = c(454,1059), breaks = seq(500,1000,100)) +
#                       theme_classic() +
#                       xlab("Wavelength (nm)") +
#                       ylab("Variance") +
#                       theme(panel.grid = element_blank(), 
#                             plot.title = element_text(face='bold', size=12),
#                             axis.text.x = element_text(face='bold', size=10, color='black'),
#                             axis.text.y = element_text(face='bold', size=10, color='black'),
#                             axis.title.y = element_text(face='bold', size=12),
#                             axis.title.x = element_text(face='bold', size=12),
#                             legend.title = element_text(face='bold', size = 10)),
#                     
#                     ggplot(SasiSDPlot) +
#                       labs(title = 'B. SWIR-Spectral') +
#                       geom_line(aes(x=Wavelength, y=Variance, group = plot_field_id, colour = PercCon)) +
#                       scale_colour_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
#                       geom_line(aes(x=Wavelength, y=Varavg), colour = 'black', linewidth = 1.1) +
#                       #annotate("rect", xmin = 899, xmax = 957, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
#                       annotate("rect", xmin = 1345, xmax = 1460, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
#                       annotate("rect", xmin = 1790, xmax = 1950, ymin = -0.001, ymax = 0.15, fill= 'grey90') +
#                       #ggbreak::scale_y_cut(breaks = 0.021, scales = 0.5, space  = .1) +
#                       scale_y_continuous(limits = c(-0.001,0.15), expand = c(0,0), breaks =seq(0.0, 0.15, 0.05)) +
#                       scale_x_continuous(limits = c(1002,2412), breaks = seq(1000,2400,400)) +
#                       theme_classic() +
#                       xlab("Wavelength (nm)") +
#                       ylab("Variance") +
#                       theme(panel.grid = element_blank(), 
#                             plot.title = element_text(face='bold', size=12),
#                             axis.text.x = element_text(face='bold', size=10, color='black'),
#                             axis.text.y = element_text(face='bold', size=10, color='black'),
#                             axis.title.y = element_text(face='bold', size=12),
#                             axis.title.x = element_text(face='bold', size=12),
#                             legend.title = element_text(face='bold', size = 10),
#                             legend.text = element_text(face = 'bold', size = 10)),
#                     
#                     nrow = 2, ncol = 1, common.legend = TRUE, legend = 'bottom')
# WavVar
# 
# ggsave(plot = WavVar, filename = 'Outputs/Figures/Wavelenght_variance.jpeg', width=3.334646, height=2, dpi=600)