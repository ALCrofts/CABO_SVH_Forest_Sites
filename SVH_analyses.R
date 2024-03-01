# Author: Anna L. Crofts, email: croftsanna@gmail.com

# Title: Code for the manuscript, "Linking aerial hyperspectral data to canopy
#        tree biodiversity: An examination of the spectral variation hypothesis."

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

# Load Data:
#   Taxonomic Dimension: species' abundance within each forest inventory plot
TaxComp <- read.csv("Data/Crofts_SVH_vegetation_survey.csv") %>%
  arrange(plot_field_id)

#   Functional Dimension: species' mean foliar trait measurements 
#FunTraits_old <- read.csv("Data/Trait_Species.csv")
FunTraits <- read.csv("~/PhD/Ch1/Trait_species_updated.csv")

#   Spectral Dimension: Continuum removed spectral reflectance within each forest inventory plot
Spectral_meta <- read.csv('Data/Crofts_SVH_spectra_metadata.csv', check.names = F)
Spectral <- read.csv("Data/Crofts_SVH_spectral_BD_dataset.csv", check.names = F)

#   - Join the additional acquisition data to the normalized spectral reflectance values.
Spectral_2 <- left_join(Spectral, # Dataframe containing spectral reflectance values.
                        Spectral_meta) %>% # Dataframe containing acquisition data. 
              select(c(
                '454.2':'2412.5', # The normalized reflectance values at each quantified wavelength band.
                spectra_BD_id, # An unique identifier for each spectral point.
                site, # The study site at which the data was acquired (i.e., either MMG or MSB).
                plot_id, # An unique identifier for each forest inventory plot. 
                plot_field_id, # Forest inventory plot id as assigned in the field.
                instrument_model # The instrument that acquired the spectral data: CASI-1500 = VNIR & SASI-640 = SWIR
              )) %>%
              relocate(c('454.2':'2412.5'), .after = instrument_model)

#   - As the two instruments were not spatially aligned, split the spectral data into two dataframes.  
Spec_VNIR <- Spectral_2 %>%
  filter(instrument_model == 'CASI-1500')

Spec_SWIR <- Spectral_2 %>%
  filter(instrument_model == 'SASI-640')
  
rm(Spectral, Spectral_2, Spectral_meta)

#   Abiotic Environment: Mean values of topographic variables for each forest inventory plot
Envr <- read.csv("Data/Crofts_SVH_envr_predictors.csv")

# Define i) plots, ii) species, and iii) sites
# i)
plot_field_id <- TaxComp$plot_field_id

# ii)
SpeciesList <-
  c(
    "Abies.balsamea..Linnaeus..Miller",
    "Acer.pensylvanicum.Linnaeus",
    "Acer.rubrum.Linnaeus",
    "Acer.saccharum.Marshall",
    "Acer.spicatum.Lamarck",
    "Alnus.incana..Linnaeus..Moench",
    "Betula.alleghaniensis.Britton",
    "Betula.papyrifera.Marshall",
    "Betula.populifolia.Marshall",
    "Carya.cordiformis..Wangenheim..K..Koch",
    "Fagus.grandifolia.Ehrhart",
    "Fraxinus.americana.Linnaeus",
    "Fraxinus.nigra.Marshall",
    "Larix.laricina..Du.Roi..K..Koch",
    "Ostrya.virginiana..Miller..K..Koch",
    "Picea.abies..Linnaeus..H..Karsten",
    "Picea.glauca..Moench..Voss",
    "Picea.rubens.Sargent",
    "Pinus.strobus.Linnaeus",
    "Populus.balsamifera.Linnaeus",
    "Populus.grandidentata.Michaux",
    "Populus.tremuloides.Michaux",
    "Prunus.pensylvanica.Linnaeus.f.",
    "Prunus.serotina.Ehrhart.var..serotina",
    "Quercus.rubra.Linnaeus",
    "Sorbus.decora..Sargent..C.K..Schneider",
    "Thuja.occidentalis.Linnaeus",
    "Tilia.americana.Linnaeus",
    "Tsuga.canadensis..Linnaeus..Carrière",
    "Ulmus.rubra.Muhlenberg"
  ) 

# iii)
SiteList <- Envr %>%
  dplyr::select(site, plot_field_id) 

# Calculate percent coniferous cover per plot
PercConifer <- TaxComp %>%
  group_by(plot_field_id) %>%
  summarise(PercCon = sum(c_across(
    cols = c(
      'Abies.balsamea..Linnaeus..Miller', # These are the 8 coniferous species observed across our study sites
      'Larix.laricina..Du.Roi..K..Koch',
      'Picea.abies..Linnaeus..H..Karsten',
      'Picea.glauca..Moench..Voss',
      'Picea.rubens.Sargent',
      'Pinus.strobus.Linnaeus',
      'Thuja.occidentalis.Linnaeus',
      'Tsuga.canadensis..Linnaeus..Carrière'
    )
  ))) %>%
  left_join(SiteList)

#write.csv(PercConifer, "Outputs/Data/PercConifer.csv", row.names = F)

# 1.0_Degree of Correspondence Analyses  ---------------------------------------------------------
#   The spectral variation hypothesis predicts that variation in plant reflectance is related to 
#   variation in plant taxonomic and functional identity (Palmer et al. 2002, Ustin & Gamon 2010).
#   Spectral and field based community properities, both composition and diversity, should be related.

## 1.1__Composition ---------------------------------------------------------
# Examine the degree of correspondence of composition using two methods:
#   i)  Multivariate correlation - ie. Procrustes analysis
#   ii) Univariate correlation   - ie. Pearson's correlation of principal component axes

### 1.1.1__Procrustes Analyses ---------------------------------------------
# Procrustes analysis is a canonical ordination method for comparing two data matrices about the same objects (here, plots). Unlike Mantel test, 
# it is appropriate to use on raw data (vs. dissimilarity matrices). The procrustes test finds a compromise ordination for two raw data matrices (by dilating, rotating, and translating),
# where the sums of squared deviations are minimized. The symmetric procrustes statistic m12^2 is a measure of similarity between the two data matrices and can be tested by the 
# Procrustean randomization test.
#     - We want to compare the PC axes that explain most of the variation in compostition (>= 90% , ie, 'total' composition)
#     - The higher value of m12^2, the weaker the relationship between the two data matrices
#     - Procrustean r = (1-SS)^1/2


#### 1.1.1.1__Explore Tax, Fun, and Spectral Comp ----------------------------
# Prior to running procrustes analyses, explore tax, fun, VNIR, and SWIR composition via examining pricipal axes of variation using PCA.

###### 1.1.1.1.1__Taxonomic Composition ----------------------------------------

# Check the structure of species abundance data frame 
sum(TaxComp == 0) # 1633 zeros
sum(TaxComp == 0) / (nrow(TaxComp) * ncol(TaxComp)) # ~78% zeros

# There is quite a lot of zeros, so apply Hellinger Transformation
#          -this transformation gives low weights to variables with low counts and many zeros
#          -square root of the relative abundance (Here, TaxComp is already characterized using relative abundances, so just need to square root transform.)
TaxComp.hellinger <- sqrt(TaxComp[, 4:33]) # Columns 4:33 contain species relataive abundances  

# Ordinate the data and examine explained variance
pca.Tax <- rda(TaxComp.hellinger) # Run PCA (rda with no explanatory variables = to pca)
summary(pca.Tax) 
#   -First two PCs explain 61.973% of the variance
#   -PC1 48.58%, PC2 13.394%, PC3 9.174% 
#   -90% variance explained by the PC1-PC8

# Visually explore ordination
plot(
  pca.Tax,
  scaling = 1,
  display = "sites",
  type = "text",
  main = "PCA for Taxonomic Composition"
)

#   Extract species scores and plot them along PCs
pca.species.scores <-
  data.frame(scores(pca.Tax, display = 'species')) %>%
  rownames_to_column('Species')

#   Visually explore main axis of tax variation
PC1.Tax <-
  ggplot(pca.species.scores, aes(x = reorder(Species, PC1), y = PC1)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Tax # Splits low and high elevation plots 

#   Visually explore second axis of tax variation
PC2.Tax <-
  ggplot(pca.species.scores, aes(x = reorder(Species, PC2), y = PC2)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Tax # Splits (kinda) Mt St Bruno and Mt Megantic


##### 1.1.1.1.2__Functional Composition ----------------------------------------
# Calculate CWM of traits for each plot
# First, select traits of interest
FunTraits <- FunTraits %>%
  dplyr::select(
    c(
      'trait_specific_leaf_area_m2_kg',
      'trait_leaf_dry_matter_content_mg_g',
      'trait_leaf_relative_water_content_perc',
      'trait_equivalent_water_thickness_cm',
      'trait_soluble_perc',
      'trait_hemicellulose_perc',
      'trait_cellulose_perc',
      'trait_lignin_perc',
      'trait_recalcitrants_perc',
      'trait_chla_mg_g_disk_mass',
      'trait_chlb_mg_g_disk_mass',
      'trait_carot_mg_g_disk_mass',
      'trait_chl_a_chl_b_ratio',
      'trait_n_perc',
      'trait_c_perc'
    )) %>%
  add_column(SpeciesList, .before = 'trait_specific_leaf_area_m2_kg')
  
# Second, calculate abundance weighted trait means
FunComp <- TaxComp %>%
  gather(species,
         abundance,
         Abies.balsamea..Linnaeus..Miller:Ulmus.rubra.Muhlenberg) %>%
  group_split(plot_field_id) %>%
  map( ~ .x %>% 
         left_join(FunTraits, by = c('species' = 'SpeciesList')), data = .x) %>%
  map_df( ~ .x %>% 
            summarise_at(vars("trait_specific_leaf_area_m2_kg":"trait_c_perc"), ~ weighted.mean(.x, abundance)
  )) %>%
  add_column(plot_field_id, .before = ('trait_specific_leaf_area_m2_kg'))

# Standardize CWMs, so traits are scaled to zero mean and unit variance
FunComp.stand <- decostand(FunComp[, 2:16], method = "standardize") 

# Ordinate the stanadardized CWMs and examine explained variance
pca.Fun <- rda(FunComp.stand)
summary(pca.Fun) 
#   - First two PCs explain 79.23% of the variance
#   - PC1 58.19%, PC2 21.04%, PC3 14.25% 
#   - >99% variance explained by the PC1-PC8

# Visually explore ordination
plot(
  pca.Fun,
  scaling = 1,
  display = "sites",
  type = "text",
  main = "PCA for Functional Composition"
)

#   Extract trait scores and plot them along PCs
pca.trait.scores <-
  data.frame(scores(pca.Fun, display = 'species')) %>%
  rownames_to_column('Traits')

#   Visually explore main axis of fun variation
PC1.Fun <-
  ggplot(pca.trait.scores, aes(x = reorder(Traits, PC1), y = PC1)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Fun # Splits conifer and deciduous chemically - cellulose perc, lignin perc, recalcitrants perc, c perc,  chla:chlb, and equivalent H2O thickness

#   Visually explore main axis of fun variation
PC2.Fun <-
  ggplot(pca.trait.scores, aes(x = reorder(Traits, PC2), y = PC2)) +
  geom_bar(stat = 'identity') +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Fun # Splits (kinda) conifer and deciduous structurally


##### 1.1.1.1.3__VNIR-Spectral Composition ----------------------------------------
# Spectral composition is defined as the mean normalized reflectance across all bands (i.e., wavelengths) per plot
VNIRComp <- Spec_VNIR %>%
  select(c(plot_field_id, '454.2':'1059.08')) %>%  # Normalized reflectance acquired btwn 454.2 nm to 1056.69 nm (or whatever wavelenghts quantified)
  group_by(plot_field_id) %>%
  summarise_all(mean)

VNIRComp <- VNIRComp[match(plot_field_id, VNIRComp$plot_field_id), ] # Make dataframe in the same order as tax and fun dataframes.
  
# Explore ordination
pca.VNIR <- rda(VNIRComp[,-1])
summary(pca.VNIR) 
#   - First two PCs explain 97.67% of the variance
#   - PC1 84.63%, PC2 13.04%, PC3 0.96% 
#   - >99% variance explained by the PC1-PC8

#   Extract wavelength scores and plot them along PCs
VNIR.band.scores <- data.frame(pca.VNIR$CA$v) %>%
  select(c(1:8)) %>% # select the first 8 PCs
  rownames_to_column(var = "Bands") %>%
  mutate(Wavelength = as.numeric(Bands), .after = Bands)

#   Visually explore main axis of VNIR variation
PC1.VNIR <- 
  ggplot(VNIR.band.scores, aes(x = Wavelength, y = PC1)) +
  geom_line() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.VNIR # Visible light, particularly blue light

#   Visually explore second axis of VNIR variation
PC2.VNIR <- 
  ggplot(VNIR.band.scores, aes(x = Wavelength, y = PC2)) +
  geom_line() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.VNIR # Blue light and NIR wavelengths

##### 1.1.1.1.3__SWIR-Spectral Composition ----------------------------------------
SWIRComp <- Spec_SWIR %>%
  select(c(plot_field_id, '972.5':'2412.5')) %>%  # Normalized reflectance acquired btwn 972.5 nm to 2412.5 nm (or whatever wavelengths were quantified)
  group_by(plot_field_id) %>%
  summarise_all(mean)

SWIRComp <- SWIRComp[match(plot_field_id, SWIRComp$plot_field_id), ] # Make dataframe in the same order as tax and fun dataframes.

# Explore ordination
pca.SWIR<-rda(SWIRComp[,-1])
summary(pca.SWIR) 
#   - First two PCs explain 86.93% of the variance
#   - PC1 58.68%, PC2 28.25%, PC3 7.23% 
#   - >99% variance explained by the PC1-PC8

#   Extract wavelength scores and plot them along PCs
SWIR.band.scores <- data.frame(pca.SWIR$CA$v) %>%
  select(c(1:8)) %>%
  rownames_to_column(var = "Bands") %>%
  mutate(Wavelength = as.numeric(Bands), .after = Bands)

#   Visually explore main axis of SWIR variation
PC1.SWIR <- 
  ggplot(SWIR.band.scores, aes(x = Wavelength, y = PC1)) +
  geom_line () +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.SWIR # Primarily longer wavelenghts > 2000 nm

#   Visually explore second axis of SWIR variation
PC2.SWIR <- 
  ggplot(SWIR.band.scores, aes(x = Wavelength, y = PC2)) +
  geom_line() +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.SWIR # Multiple peaks across SWIR

# # Remove exploratory objects not needed going forward
# rm(PC1.Tax, PC2.Tax, PC1.Fun, PC2.Fun, PC1.VNIR, PC2.VNIR, PC1.SWIR, PC2.SWIR)

### 1.1.1.2_Run Procrustes Analyses -----------------------------------------

# PCs 1-8 contains greater than 90% of variation in composition across all dimensions 
#   Create dataframes of site scores across PCs 1-8.

# PCA_sitescores function ordinates data and creates a dataframe with PC1-8
PCA_sitescores <- function(Comp, # Composition objects (site by species / CWM traits / avg. reflectance)
                           plot_field_id, # Vector of plot_ids
                           dimension) { # Type of composition data (i.e., Taxonomic, Functional, Spectral)
  PCA <- rda(Comp)
  PCA_outputs <- data.frame(PCA$CA$u) %>%
    select(c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8')) %>%
    cbind(plot_field_id) %>%
    mutate(dimension = dimension)
  
  print(PCA_outputs)
}

# Create a dataframes with site scores of PC1-8 for each dimension
TaxPCA_scores <- PCA_sitescores(TaxComp.hellinger, plot_field_id, 'Taxonomic')
FunPCA_scores <- PCA_sitescores(FunComp.stand, plot_field_id, 'Functional')
VNIRPCA_scores <- PCA_sitescores(VNIRComp[,-1], plot_field_id, 'VNIR')
SWIRPCA_scores <- PCA_sitescores(SWIRComp[,-1], plot_field_id, 'SWIR')

PCA_outputs <- rbind(TaxPCA_scores, FunPCA_scores, VNIRPCA_scores, SWIRPCA_scores)
write.csv(PCA_outputs, "Outputs/Data/PCA_outputs_updated.csv", row.names = F)

# Run procrustes analyses and procrustean permutation test (n = 999)
Tax.VNIR.out <- # Output of procustes analysis btwn taxonomic and VNIR-spectral composition
  protest(
    X = TaxPCA_scores[, 1:8], # Taxonomic site scores for PCs 1-8
    Y = VNIRPCA_scores[, 1:8], # VNIR-spectral site scores for PCs 1-8
    symmetric = T, 
    scores = "sites",
    permutations = 999
  )

Tax.SWIR.out <- # Output of procustes analysis btwn taxonomic and SWIR-spectral composition
  protest(
    X = TaxPCA_scores[, 1:8], # Taxonomic site scores for PCs 1-8
    Y = SWIRPCA_scores[, 1:8], # SWIR spectral site scores for PCs 1-8
    symmetric = T,
    scores = "sites",
    permutations = 999
  )

Fun.VNIR.out <- # Output of procustes analysis btwn functional and VNIR-spectral composition
  protest(
    X = FunPCA_scores[, 1:8], # Functional site scores for PCs 1-8
    Y = VNIRPCA_scores[, 1:8], # VNIR-spectral site scores for PCs 1-8
    symmetric = T,
    scores = "sites",
    permutations = 999
  )

Fun.SWIR.out <- # Output of procustes analysis btwn functional and SWIR-spectral composition
  protest(
    X = FunPCA_scores[, 1:8], # Functional site scores for PCs 1-8
    Y = SWIRPCA_scores[, 1:8],  # SWIR spectral site scores for PCs 1-8
    symmetric = T,
    scores = "sites",
    permutations = 999
  )

# Create a dataframe detailing procustes output
Procrustes <- data.frame(Field = c('Taxonomic', 'Taxonomic', 'Functional', 'Functional'),
                         Spectral = c('VNIR', 'SWIR', 'VNIR', 'SWIR'),
                         r =c(Tax.VNIR.out$t0, Tax.SWIR.out$t0, Fun.VNIR.out$t0, Fun.SWIR.out$t0),
                         m12.squared =c(Tax.VNIR.out[["ss"]], Tax.SWIR.out[["ss"]], Fun.VNIR.out[["ss"]], Fun.SWIR.out[["ss"]]),
                         p = c(Tax.VNIR.out[["signif"]], Tax.SWIR.out[["signif"]], Fun.VNIR.out[["signif"]], Fun.SWIR.out[["signif"]]))

write.csv(Procrustes, "Outputs/Statistics/Procrustes_outputs.csv")

# Plot fit of procrustes analyses  
ProcrustesFit <-
  ggplot(Procrustes, aes(x = factor(Spectral, levels = c('VNIR', 'SWIR')), 
                         y = r, 
                         fill = factor(Field, levels =c("Taxonomic", "Functional")))) +
  geom_bar(stat = 'identity', position = position_dodge(), colour = 'black') +
  scale_fill_manual(values = c("Grey40", "Grey80")) +
  geom_text(aes(label = round(r,3)), position=position_dodge(width=0.9), vjust=-0.25, size=4)+
  scale_y_continuous(expand=c(0,0),breaks=seq(0, 1, by=0.25), limits=c(0,1))+
  theme_classic() +
  ylab("Procrustean r")+
  xlab("") +
  theme(axis.text.x=element_text(face='bold', size=12, colour = 'black'),
        axis.text.y=element_text(face='bold', size=10, colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.title.y =element_text(face='bold', size=12),
        legend.text = element_text(face='bold', size=12),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.95))
ProcrustesFit
ggsave("Outputs/Figures/ProcrustesFit.jpeg", width=6, height=4, dpi=300)


# rm(Tax.VNIR.out, Tax.SWIR.out, Fun.VNIR.out, Fun.SWIR.out, ProcrustesFit, Procrustes)

### 1.1.2_Pearson's Correlations ----------------------------------------------------
#   Examine if the individual PCs are associated with each other.
#   Visualize species, trait, wavelength scores along PC axes


#### 1.1.2.1__Run PC Correlations --------------------------------------------

# PCA_outputs <- read.csv("Outputs/Data/PCA_outputs.csv")

PCA_outputs2 <- PCA_outputs %>%
  pivot_wider(names_from = dimension,
              values_from = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8))

# Pearson's correlations of the first 3 PC axes
PCACor.Total <- cor(PCA_outputs2[, 2:13]) # Columns 2-13 are characterized by the first 3 PC axes for each dimension (i.e., tax, fun, VNIR-comp, SWIR-comp) 

# Name the rows and columns of correlation matrix
colnames(PCACor.Total) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                            "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                            "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")
rownames(PCACor.Total) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                            "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                            "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")

# Significance text to produce p-values for each pair of correlations 
PCACor.Test <- cor.mtest(PCA_outputs2[,2:13], conf.level = 0.95)

# Name the rows and columns of the p-value matrix
colnames(PCACor.Test[['p']]) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                                  "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                                  "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")
rownames(PCACor.Test[['p']]) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                                  "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                                  "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3")

# Create a dataframe with correlation coefficients for the first 3 PCs across dimensions
PCACor.Total2 <- pull_lower_triangle(PCACor.Total) %>%
  gather(Var2, r, 2:13) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) %>%
  mutate(r = round(r, 2)) 

# Create a dataframe with correlation p-values for the first 3 PCs across dimensions
PCACor.Test2 <- pull_lower_triangle(PCACor.Test[["p"]]) %>% 
  gather(Var2, p, 2:13) %>%
  rename(Var1 = rowname) %>%
  mutate(p = as.numeric(p)) %>%
  filter(!is.na(p)) %>%
  mutate(p = round(p, 3)) 

# Merge correlation coefficients and p-values and filter to only keep field-based composition to spectral composition comparisions
PCACor.Total3 <- left_join(PCACor.Total2, PCACor.Test2) %>%
  filter(Var1 == "VNIR_PC1" | Var1 == "VNIR_PC2" | Var1 == "VNIR_PC3" |
         Var1 == "SWIR_PC1" | Var1 == "SWIR_PC2" | Var1 == "SWIR_PC3" ) %>%
  filter(Var2 == "Tax_PC1" | Var2 == "Tax_PC2" | Var2 == "Tax_PC3" |
         Var2 == "Fun_PC1" | Var2 == "Fun_PC2" | Var2 == "Fun_PC3" ) %>%
  rename(Spectra = Var1, Field = Var2)

write.csv(PCACor.Total3, "Outputs/Statistics/PCA_correlations_outputs.csv")

# Plot the Correlation matrix, Figure S5  
PCACor.Total4 <- left_join(PCACor.Total2, PCACor.Test2) %>%
  mutate(r2 = ifelse(p <= 0.05, r, NA)) %>%  # Remove r values for non-signif correlations -- FOR PLOTTING PURPOSES ONLY!!
  rename(Spectra = Var1, Field = Var2)

PCACorTotalMatrix <- 
  ggplot(data = PCACor.Total4, aes(Field, 
                                   Spectra, 
                                   fill = r))+
  scale_x_discrete(limits = c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1",
                              "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2",
                              "Tax_PC3", "Fun_PC3", "VNIR_PC3")) +
  scale_y_discrete(limits = c("SWIR_PC3","VNIR_PC3","Fun_PC3", "Tax_PC3",
                              "SWIR_PC2","VNIR_PC2","Fun_PC2", "Tax_PC2",
                              "SWIR_PC1","VNIR_PC1","Fun_PC1")) +
  geom_tile(color = "white")+
  labs(title = "") +
  geom_text(aes(Field, Spectra, label = r), color = "black", size = 4) +
  geom_text(aes(Field, Spectra, label = r2), na.rm = T, colour = 'black', fontface = 'bold', size = 4) +
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

ggsave(plot = PCACorTotalMatrix, filename = 'Outputs/Figures/FigS5_PCACor_Total.jpeg', width = 6, height = 6, dpi = 300)

# For presentations, make a graph that depicts the PC correlation fit in an analogous way to the procrustes figure.
# PCs <- ggplot(PCACor.Total4 %>% filter(Var1 == "VNIR_PC1" & Var2 == "Tax_PC1" |
#                                        Var1 == 'VNIR_PC1' & Var2 == 'Fun_PC1' ),
#               aes(x = factor(Var2, levels = c("Tax_PC1", "Fun_PC1"), labels = c("Taxonomic", 'Functional')),
#                   fill = factor(Var2, levels = c("Tax_PC1", "Fun_PC1"), labels = c("Taxonomic", 'Functional')),
#                   y = r)) +
#   geom_bar(stat = 'identity', position = position_dodge(), colour = 'black') +
#   scale_fill_manual(values = c("Grey40", "Grey80")) +
#   geom_text(aes(label = r), position=position_dodge(width=0.9), vjust=-0.25, size=4)+
#   scale_y_continuous(expand=c(0,0),breaks=seq(0, 1, by=0.25), limits=c(0,1))+ 
#   theme_classic() +
#   ylab("Correlation with Spectral Composition")+
#   xlab("")+
#   theme(panel.grid = element_blank(),
#         axis.text.x=element_text(face='bold', size=12, colour = 'black'),
#         axis.ticks.x = element_blank(),
#         axis.text.y=element_text(face='bold', size=10),
#         axis.title.y =element_text(face='bold', size=12),
#         legend.text = element_text(face='bold', size=12),
#         legend.title = element_blank(),
#         legend.position = 'none') 
# PCs
# ggsave("Outputs/Figures/PCs_pres.jpeg", width=5, height=4.5, dpi=600)

rm(PCACor.Total, PCACor.Total2, PCACor.Total3, PCACor.Total4, PCACor.Test, PCACor.Test2, PCACorTotalMatrix)

#### 1.1.2.2__Understand PC axes -------------------------------------
#   - How do the primary principal components relate to the temperate-to-boreal gradient? 
#   - How are species, traits, and wavelengths distributed along the primary principal components?

# Pearson's correlations of primary PC axes and percent conifer content
PCA_outputs3 <- PCA_outputs2 %>%
  left_join(PercConifer)

#  Correlate PC1 of each dimension (cols 2:5) and percent coniferous cover (column 34)
PCACor.PC1 <- cor(PCA_outputs3[, c(2:5,34)]) 
colnames(PCACor.PC1) <- c("Taxonomic", "Functional", "VNIR-Spectral", "SWIR-Spectral", "Perc. Coniferous")
rownames(PCACor.PC1) <- c("Taxonomic", "Functional", "VNIR-Spectral", "SWIR-Spectral", "Perc. Coniferous")

PCACor.PC1_2 <- pull_lower_triangle(PCACor.PC1) %>% 
  gather(Var2, r, 2:6) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) %>%
  mutate(r = round(r, 2)) 

#  Visualize associations 
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
  summarise_at(.vars = c(4:33), # Columns 4:33 contain species abundances
               .funs = sum) %>%
  gather(key = 'Species', value = 'Abundance') %>%
  arrange(desc(Abundance)) %>%
  slice(1:15) # Only retain the 15 most abundant species 

pca.species.scores2 <- pca.species.scores %>%
  filter(Species %in% TopSpecies$Species) %>% # Only keep 15 most abundant species to simplify figure
  mutate(Names = c('Abies balsamea', 'Acer pensilvanicum', 'Acer rubrum', 'Acer saccharum', 'Betula alleganiensis', 'Betula papyrifera', 'Betula popufolia', 'Fagus grandifolia','Ostrya virginiana','Picea rubens', 'Populus tremuloides','Prunus strobus','Quercus rubra','Tilia americana','Tsuga canadensis'))

PC1.Tax <- ggplot(pca.species.scores2, aes(x=reorder(Names, PC1), y=PC1)) +
  geom_bar(stat='identity', colour = 'black', fill = 'grey60') +
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
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color ='black', face = 'bold'),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 12),
        legend.title = element_text(face = 'bold', size = 10)) 
PC1.Tax #Splits low and high elevation plots 

pca.trait.scores2 <- pca.trait.scores %>%
  mutate(Traits2 = c('SLA', 'LDMC','RWC', 'EWT', 'Soluable', 'Hemicellulose', 'Cellulose', 'Lignin', 'Recalcitrants', 'Chl a', 'Chl b', 'Carotenoids', 'Chl a : Chl b', 'N', 'C'))

PC1.Fun <- ggplot(pca.trait.scores2, aes(x=reorder(Traits2, PC1), y=PC1)) +
  geom_bar(stat='identity', colour = 'black', fill = 'grey60') +
  xlab("")+
#  ylab("Functional PC1 \n (58%)") +
  ylab("PC1 Scores") +
  scale_x_discrete(labels = c("SLA", "Carot.", "Chl b", "Chl a", "RWC", "N", "Hemi.", "LDMC", "Sol.", "Cell.", "Chla:Chlb", "Lignin","Recalc.", "C", "EWT")) +
  labs(title = 'Functional', subtitle = 'Variance Explained = 58%') +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color='black', face = 'bold'),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        legend.position = 'none') 
PC1.Fun #Splits conifer and deciduous chemically - cellulose perc, lignin perc, recalcitrants perc, c perc,  chla:chlb, and equivalent H2O thickness  

VNIR.band.scores2 <- VNIR.band.scores %>%
  mutate(Region = ifelse(Wavelength <= 500, 'blue',
                         ifelse(Wavelength > 500 & Wavelength <= 600, 'green',
                                ifelse(Wavelength > 600 & Wavelength <= 730, 'red',
                                       ifelse(Wavelength > 730 & Wavelength <=780, 'rededge', 'NIR'))))) %>%
  full_join(data.frame(Region = 'SWIR',
                       PC1 = NA,
                       Wavelength = 454))

PC1.VNIR <- ggplot(VNIR.band.scores2, aes(x = Wavelength, y = PC1)) +
  geom_area(fill = 'steelblue4') +
  geom_area(mapping = aes(x = ifelse(Wavelength > 590, Wavelength, 0)), fill = 'mediumseagreen') +
  geom_area(mapping = aes(x = ifelse(Wavelength > 700, Wavelength, 0)), fill = 'tomato3') +
  geom_area(aes(fill = Region)) +
  geom_line() +
  scale_fill_manual(values = c('steelblue4', 'mediumseagreen','tomato3','lightcoral',  'indianred3',  'rosybrown3'),
                    breaks = c('blue','green','red','rededge','NIR','SWIR'),
                    name = 'Spectral Region', 
                    labels = c('Blue','Green','Red', 'Red Edge', 'NIR', 'SWIR'),
                    limits = c('blue','green','red','rededge','NIR','SWIR'),
                    guide = guide_legend(title.hjust = 0.5, title.position = "top", nrow = 1)) + 

  annotate('rect', xmin = 899, xmax = 957, ymin = -0.21, ymax = 0.21, colour = 'grey90', fill = 'grey90') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab("Wavelength (nm)") + 
#  ylab("VNIR-Spectral PC1 \n (85%)") + 
  ylab('PC1 Scores') +
  labs(title = 'VNIR-Spectral', subtitle = 'Variance Explained = 84%') +
  scale_y_continuous(limits = c(-0.21, 0.21), expand = c(0, 0)) +
  scale_x_continuous(limits = c(454, 1059), breaks = seq(500, 1000, 100)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color ='black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5)
PC1.VNIR 

SWIR.band.scores2 <- SWIR.band.scores %>%
  mutate(Region = ifelse(Wavelength <= 1400, 'NIR', 'SWIR'))

PC1.SWIR <- ggplot(SWIR.band.scores2, aes(x = Wavelength, y = PC1, fill = Region)) +
  geom_area() +
  geom_line() +
  scale_fill_manual(values = c('indianred3', 'rosybrown3')) +
  annotate('rect', xmin=c(1333,1790), xmax=c(1466,1950), ymin=c(-0.33,-0.33), ymax=c(0.33,0.33), colour = 'grey90', fill = 'grey90') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab("Wavelength (nm)") +
  #  ylab("SWIR-Spectral PC2 \n (26%)") +
  ylab("PC2 Scores") +
  labs(title = 'SWIR-Spectral', subtitle = 'Variance Explained = 59%') +
  scale_y_continuous(limits = c(-0.33, 0.33), breaks = seq(-0.3, 0.3, 0.15), expand = c(0,0)) +
  scale_x_continuous(limits = c(1002,2412), breaks = seq(1000,2400,200)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color ='black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5)
PC1.SWIR

PC2.SWIR <- ggplot(SWIR.band.scores2, aes(x=Wavelength, y=PC2, fill = Region)) +
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
        legend.text = element_text(size = 10),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color ='black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5)
PC2.SWIR


windows(width = 7.08661, height = 6)
PCPlots <- ggarrange(PC1.Tax, PC1.Fun, PC1.VNIR, PC2.SWIR,
          nrow = 2, ncol = 2,
          common.legend = T, legend = 'bottom',
          labels =  'AUTO')
PCPlots
ggsave(plot = PCPlots, filename = 'Outputs/Figures/Fig2_PC_plots.pdf', width = 7.08661, height = 6, dpi = 600)


# rm(PCA_outputs, PCA_outputs2, PCA_outputs3, PCACor.PC1, PCACor.PC1_2, PC1CorMatrix,
#    TaxComp.hellinger, TaxPCA_scores, pca.Tax, pca.species.scores,  pca.species.scores2, TopSpecies,
#    FunComp.stand, FunPCA_scores,pca.Fun, pca.trait.scores, pca.trait.scores2,
#    VNIRPCA_scores, pca.VNIR, VNIR.band.scores, VNIR.band.scores2,
#    SWIRPCA_scores, pca.SWIR, SWIR.band.scores, SWIR.band.scores2,
#    PC1.Tax, PC1.Fun, PC1.SWIR, PC2.Tax, PC2.Fun, PC2.VNIR, PC2.SWIR, PCPlots)

# 1.2__Diversity  --------------------------------------------------------------
#   - First, calculate diversity indices for each biodiversity dimension
#   - Then, examine the degree of association between spectral and taxonomic/functinal indices
#       * 1st via linear correlation for all indices 
#       * 2nd via non-linear association for subset (?) - best correlated indices


## 1.2.1__Calculate Diversity Indices --------------------------------------


### 1.2.1.1__Taxonomic Diversity Indices  -------------------------------------------------------------------------
# Calculate hill diversity for each field plot:
#   q0 = species richness, the number of species per plot
#   q1 = Shannon diversity, the uncertainty in predicting the identity of a new indiv. given species richness and abundance
#   q2 = Simpson diversity, the prob. that two indiv. drawn at random will be different species
#   J  = Pielou's eveness, how close in abundance each species is

TaxDiv <- TaxComp %>%
  summarise(plot_field_id = plot_field_id,
            q0 = hill_taxa(TaxComp[, 4:33], q = 0, MARGIN = 1),  # Cols 4:33 contain species abundances
            q1 = hill_taxa(TaxComp[, 4:33], q = 1, MARGIN = 1),  # Cols 4:33 contain species abundances
            q2 = hill_taxa(TaxComp[, 4:33], q = 2, MARGIN = 1)) %>% # Cols 4:33 contain species abundances
  mutate(J =  log(q1)/log(q0)) %>%
  mutate(q1 = ifelse(plot_field_id == 'Plot13_Maples', NA, q1), # Plot13_Maples only has 1 species, can't calculate exp. Shannon's index
         q2 = ifelse(plot_field_id == 'Plot13_Maples', NA, q2), # Plot13_Maples only has 1 species, can't calculate inverse Simpson's index
         J  = ifelse(plot_field_id == 'Plot13_Maples', NA, J))  # Plot13_Maples only has 1 species, can't calculate Pielou's evenness index

### 1.2.1.2__Functional Diversity Indices ----------------------------------------
# Use FD package to calculate functional diversity indices for each plot:
#   FRic = Functional Richness, convex hull volume of functional trait space
#   FEve = Functional Evenness, quantifies the regularity of interspecies distances and the homogeneity of species abundances (btwn 0-1)
#   FDiv = Functional Divergence, the deviation of species from the centre of gravity functional trait space
#   FDis = Functional Dispersion, the mean distance of a species from the centroid of all species in functional trait space
#   RaoQ = Rao's Quadratic Entropy, the sum of distances of btwn species pairs wighted by their relative abundance 


dbFD <- dbFD(
            FunTraits %>% 
              column_to_rownames(var = 'SpeciesList'), # Species-mean trait values dataframe
            TaxComp[, 4:33], # Site-by-species dataframe; Cols 4:33 contain species abundances
            w.abun = T, 
            stand.x = T, 
            calc.FRic = T,
            calc.FDiv = T
          )

FunDiv <- data.frame(plot_field_id = plot_field_id,
                     FRic = dbFD[["FRic"]],
                     FEve = dbFD[["FEve"]],
                     FDiv = dbFD[['FDiv']],
                     FDis = dbFD[['FDis']],
                     RaoQ = dbFD[['RaoQ']])

### 1.2.1.3__VNIR-Spectral Diversity Indices ----------------------------------------------
# Calculate three spectral diversity metrics for each plot: 
#   CV  = coefficient of variation, the average coefficient of variation for all wavelengths of the data within each plot. Larger CV corresponds to higher SD.
#   CHV = convex hull volume, the volume of a polygon of pixels/points forming a convex hull of the first three PC. Larger CHV corresponds with higher SD.
#   SV = spectral variance, the total sums of squared deviation for a community standardized by number of spectral points in the community (Laliberte et al. 2020). Higher SDa (aka neighboring pixels are spectral different from one another) corresponds with greater diversity. 

# calc.CV calculates the spectral diversity metric coefficient of variaiton (CV)
calc.CV <- function(spectral_df, # Dataframe with spectra reflectance values 
                    areas_of_interest, # What you want to calculate spectral diversity for, here it's plots.
                    wavelengths, # Cols where spectral reflectance values are
                    rarefraction, # If TRUE, spectral observations are standardized and randomly resampled. If FALSE, uses all spectral observations as is. 
                    n # Number of random resampling events, if rarefraction = T.
                    ){ 
  
  if(rarefraction == TRUE){
    
  # a) Determine min number of spectral points/pixels observed in a plot
  spectral_points <- spectral_df %>% 
    group_by({{areas_of_interest}}) %>%
    summarise(points = n())
  
  min_points = min(spectral_points$points)
  
  # b) randomly resample plots to minimum number of observed spectral points/pixels
  spectral_df2 <- replicate(n, # replicate dataframe n times
    spectral_df %>%
    select(c({{areas_of_interest}}, {{wavelengths}})) %>%
    group_by({{areas_of_interest}}) %>%
    sample_n(min_points, replace = F), # sample n points, where n = min_points, without replacement
    F) # FALSE gives a list output 
  
  # c) calculate CV for each plot for each resampling event
  CV <- spectral_df2 %>%
    map(~ .x %>% group_by({{areas_of_interest}}) %>%     
          summarise_all(~sd(.x)/abs(mean(.x))) %>% # calculates CV for each wavelength
          rowwise({{areas_of_interest}}) %>%
          summarise(CV = sum(c_across(cols = everything()), na.rm = T) / (ncol(.) - sum(is.na(c_across(everything())))))) %>% # sums across wavelengths and standardizes by the number of bands
    do.call(rbind, .) %>% # collapse list of resampling events into a dataframe
    group_by({{areas_of_interest}}) %>% # group resampling events by plot
    summarise(CV = mean(CV)) # calculate average CV for each plot across the resampling events

  return(CV)
  }
  
  if(rarefraction == FALSE){
    
    CV2 <- spectral_df %>%
      select(c({{areas_of_interest}}, {{wavelengths}})) %>%
      group_by({{areas_of_interest}}) %>%
      summarise_all(~sd(.)/abs(mean(.))) %>%
      rowwise({{areas_of_interest}}) %>%
      summarise(CV = sum(c_across(cols = everything()), na.rm = T) / (ncol(.) - sum(is.na(c_across(everything())))))
    
    return(CV2)
  }
}

# calc.CHV calculates the spectral diversity metric convex hull volume (CHV)
calc.CHV <- function(spectral_df, # Dataframe with spectral reflectance values
                     areas_of_interest, # What you want to calculate spectral diversity for, ie. grouping variable, here it is plots. 
                     wavelengths, # Cols where spectral reflectance values are
                     rarefraction, # If TRUE, spectral observations are standardized and randomly resampled. If FALSE, uses all spectral observations as is. 
                     n # Number of random resampling events, if rarefraction = T.
                     ){  
  
  if(rarefraction == TRUE){
    
  # a) Determine min number of spectral points/pixels observed in a plot  
  spectral_points <- spectral_df %>% 
    group_by({{areas_of_interest}}) %>%
    summarise(points = n())
  
  min_points = min(spectral_points$points)
  
  # 1b) Ordinate the data  
  PCA <- spectral_df %>%
    select(c({{wavelengths}})) %>%
    rda(scale = F) 
  
  # 2b) Create a dataframe of the first 3 PC axes
  PCA_scores <- data.frame(PCA$CA$u) %>%
    select(c(1:3)) %>%
    cbind({{spectral_df}}) %>%
    select(c('PC1', 'PC2', 'PC3', {{areas_of_interest}}))
  
  # c) Convert to quoted areas_of_interest object, needed to name list objects created in next step.
   areas_of_interest2 <- deparse(substitute(areas_of_interest))
   
  # d) calculate average CHV for each plot across all resampling events
  CHV <- PCA_scores %>%
    group_split({{areas_of_interest}}) %>% # split by plot
    set_names(map(., ~unique(.[[areas_of_interest2]]))) %>% # assign plot names to list objects 
    map(~replicate(n,
                   .x %>%
                   group_by({{areas_of_interest}}) %>%
                   sample_n(min_points, replace = F), F) %>% # resample PC scores n-times and retain the minimum number of observed spectral observations, where resampling occurs without replacement  
        map(~convhulln(.x[-4], option = 'FA')) %>% # fit convex hulls to each plot and each resampling event
        map_dbl('vol')) %>% # retain the convex hull volume for each plot and each resampling event
    map_df(~mean(.x)) %>% # calculate mean CHV for each plot across all resampling events
    do.call(rbind, .) %>% # collapse to dataframe of plots and CHV
    as.data.frame() %>%
    rownames_to_column('plot_field_id') %>%
    rename(CHV = V1)
  
  rm(areas_of_interest2)
  
  return(CHV)
  }
  
  if(rarefraction == FALSE){
  
  areas_of_interest2 <- deparse(substitute(areas_of_interest))
    
  PCA2 <- spectral_df %>%
    select(c({{wavelengths}})) %>%
    rda(scale = F) 
  
  CHV2 <- data.frame(PCA2$CA$u) %>%
    select(c(1:3)) %>%
    cbind({{spectral_df}}) %>%
    select(c('PC1', 'PC2', 'PC3', {{areas_of_interest}})) %>%
    group_split({{areas_of_interest}}) %>%
    set_names(map(., ~unique(.[[areas_of_interest2]]))) %>%
    map(~convhulln(.x[-4], option = 'FA')) %>%
    map_dbl('vol') %>%
    as.data.frame() %>%
    rownames_to_column('plot_field_id') %>%
    rename(CHV = '.')
  
  rm(areas_of_interest2) 
  
  return(CHV2)
  }
  
} 

# calc.SV calculates the spectral diversity metric spectral variation (SV)
calc.SV <- function(spectral_df, # Dataframe with spectra reflectance values 
                    areas_of_interest, # What you want to calculate spectral diversity for, here it's plots.
                    wavelengths, # Cols where spectral reflectance values are
                    rarefraction,  # If TRUE, spectral observations are standardized and randomly resampled. If FALSE, uses all spectral observations as is. 
                    n # Number of random resampling events, if rarefraction = T.
                    ){
  
  # a) Determine the number of spectral observations per plot
  spectral_points <- spectral_df %>% 
    group_by({{areas_of_interest}}) %>%
    summarise(points = n())
  
  if(rarefraction == TRUE){
    
    # b) Determine the minimum number of spectral observations observed
    min_points = min(spectral_points$points) 
    
    # c) randomly resample plots to minimum number of observed spectral points/pixels
    spectral_df2 <- replicate(n, # replicate resampling n-times
        spectral_df %>%
        select(c({{areas_of_interest}}, {{wavelengths}})) %>%
        group_by({{areas_of_interest}}) %>%
        sample_n(min_points, replace = F), # sample n points, where n = min_points, without replacement
        F) # FALSE returns a list
    
    SV <- spectral_df2 %>%
      map(~.x %>%
            group_by({{areas_of_interest}}) %>%
            summarise_all(~sum((.x - mean(.x))^2)) %>% # calculate the sum of the squared difference from the mean for each wavelength
            rowwise({{areas_of_interest}}) %>%
            summarise(SS = sum(c_across(cols = everything()))) %>% # calculate the sum of squares across all wavelengths
            summarise(SV = SS / (min_points - 1))) %>% # divide the sum of squares by the number of spectral observations (here, it equals the min points observed)
      do.call(rbind, .) %>%
      group_by({{areas_of_interest}}) %>%
      summarise(SV = mean(SV))
    
    return(SV)
  }
  
  if(rarefraction == FALSE){
    
    SV2 <- spectral_df %>%
      select(c({{wavelengths}}, {{areas_of_interest}})) %>%
      group_by({{areas_of_interest}}) %>%
      summarise_all(~sum((.x - mean(.x))^2)) %>%
      rowwise({{areas_of_interest}}) %>%
      summarise(SS = sum(c_across(cols = everything()))) %>%
      left_join(spectral_points) %>%
      summarise(SV = SS / (points - 1)) # divide by the number of spectral observations, here it can vary by plot
    
    return(SV2)
  }
} 

# Calculate spectral diversity metrics for VNIR-spectral
VNIRDiv <- calc.CV(Spec_VNIR, plot_field_id, c("454.2": "1059.08"), TRUE, 999) %>%
  left_join(calc.CHV(Spec_VNIR, plot_field_id, c("454.2": "1059.08"), TRUE, 999)) %>%
  left_join(calc.SV(Spec_VNIR, plot_field_id, c("454.2": "1059.08"), FALSE, NA)) %>%
  rename_at(.vars = c('CV', 'CHV', 'SV'), .funs = ~paste0(.x, '_VNIR'))

### 1.2.1.4__SWIR-Spectral Diversity Indices ----------------------------------------------
# Calculate spectral diversity metrics for VNIR-spectral
SWIRDiv <- calc.CV(Spec_SWIR, plot_field_id, c("972.5": "2412.5"), TRUE, 999) %>%
  left_join(calc.CHV(Spec_SWIR, plot_field_id, c("972.5": "2412.5"), TRUE, 999)) %>%
  left_join(calc.SV(Spec_SWIR, plot_field_id, c("972.5": "2412.5"), FALSE, NA)) %>%
  rename_at(.vars = c('CV', 'CHV', 'SV'), .funs = ~paste0(.x, '_SWIR'))

## 1.2.2__Diversity Correlations ----------------------------------------------
#   Examine the association of between spectral diversity indices and taxonomic / function diversity indices.
#     Pearson correlation measures the strength of the linear relationship between two variables

# First, create a dataframe containing all diversity indices
DivIndices <- TaxDiv %>%
  left_join(FunDiv) %>%
  left_join(VNIRDiv) %>%
  left_join(SWIRDiv) 

write.csv(DivIndices, file = "Outputs/Data/DiversityIndices.csv", row.names = F)

rm(TaxDiv, dbFD, FunDiv, VNIRDiv, SWIRDiv)

# Do correlations and keep correlation statistics between spectral div indices and tax/fun div indices
#   - Pearson's correlation
#   - Complete.obs removes NA values by case wise deletion

# Do correlation, where we exclude Na values (ie. plots where div index could be calculated)
DivCor.total <- cor(DivIndices[,-1], use="complete.obs")

# Name matrix cols and rows
colnames(DivCor.total) <- c("q0", "q1", "q2", "J",
                            "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                            "CV_VNIR", "CHV_VNIR", "SV_VNIR",
                            "CV_SWIR", "CHV_SWIR", "SV_SWIR")

rownames(DivCor.total) <- c("q0", "q1", "q2", "J",
                            "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                            "CV_VNIR", "CHV_VNIR", "SV_VNIR",
                            "CV_SWIR", "CHV_SWIR", "SV_SWIR")

# Test significance of correlations
DivCor.Test = cor.mtest(DivIndices[,-1], use="complete.obs", conf.level = 0.95)

colnames(DivCor.Test[['p']]) <- c("q0", "q1", "q2", "J",
                                  "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                                  "CV_VNIR", "CHV_VNIR", "SV_VNIR",
                                  "CV_SWIR", "CHV_SWIR", "SV_SWIR")

rownames(DivCor.Test[['p']]) <- c("q0", "q1", "q2", "J",
                                  "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                                  "CV_VNIR", "CHV_VNIR", "SV_VNIR",
                                  "CV_SWIR", "CHV_SWIR", "SV_SWIR")

# Create dataframe of pearson's correlation coeficients
DivCor.Total2 <- pull_lower_triangle(DivCor.total) %>% 
  gather(Var2, r, 2:16) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) 

# Create dataframe of significance values
DivCor.Test2 <- pull_lower_triangle(DivCor.Test[["p"]]) %>% 
  gather(Var2, p, 2:16) %>%
  rename(Var1 = rowname) %>%
  mutate(p = as.numeric(p)) %>%
  filter(!is.na(p))  

# Only keep spectral - field relationships
DivCor <- left_join(DivCor.Total2, DivCor.Test2) %>%
  filter(Var1 == "CV_VNIR" | Var1 == "CHV_VNIR" |Var1 == "SV_VNIR" |
         Var1 == "CV_SWIR" | Var1 == "CHV_SWIR" |Var1 == "SV_SWIR") %>%
  filter(Var2 == "q0" | Var2 == "q1" | Var2 == "q2" | Var2 == "J" |
         Var2 == "FRic" | Var2 == "FEve" | Var2 == "FDiv" | Var2 == "FDis" | Var2 == "RaoQ") %>%
  rename(Spectral = Var1, Field = Var2)

write.csv(DivCor, file = "Outputs/Statistics/DiversityCorrelations.csv", row.names = F)

# Summarize correlations by Wavelength Region, Field Dimension, and Spectral Div Metric
DivCor2 <- DivCor %>%
  separate(Spectral, c("Spec_metric","Spec_region"), sep = '[_]') %>%
  mutate(Field_dimension = ifelse(Field == "q0" |
                                  Field == "q1" | 
                                  Field == "q2" | 
                                  Field == "J", "Tax", "Func"),
         .after = Field) %>%
  rename(Field_metric = Field)

# Calculate mean (and se) correlation coefficient for spectral metrics 
DivCor2_summary <- DivCor2 %>%
  group_by(Spec_region, Field_dimension, Spec_metric) %>%
  mutate(mean = mean(r), se = sd(r)/sqrt(n()))

# Create wide format, Table 2
DivCor3 <- DivCor2 %>%
  mutate(r = round(r, 3)) %>%
  mutate(p = round(p, 3)) %>%
  pivot_wider(names_from = Spec_metric, values_from = c("r","p")) 

write.csv(DivCor3, "Outputs/Statistics/Table2_DiversityCorrelations.csv")

# Plot mean correlation with spectral diversity, Figure 3 
windows(width = 3.34646, height = 3.34646)

DivCor_plot <- ggplot(DivCor2_summary, aes(x = factor(Field_dimension, level = c('Tax','Func')), 
                                           y = mean, 
                                           group = factor(Spec_metric, level = c('CV','CHV',"SV")))) +
  geom_bar(stat = 'identity', position = position_dodge(), 
           aes(fill = factor(Spec_metric, level = c('CV', 'CHV', 'SV'))), colour = 'black') +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.9), width = 0.2) +
  facet_wrap(~factor(Spec_region, level = c('VNIR','SWIR'))) +
  scale_fill_grey(start = 0.8, end = 0.2, name = 'Spectral Metric') +
  scale_x_discrete("", labels = c("Taxonomic", 'Functional')) +
  scale_y_continuous("Mean Correlation with Spectral Diversity \n(Pearson's r)", expand = c(0,0), limits = c(-0,1), breaks = c(0, .25, .5, .75, 1))+
  theme_classic() +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(face = 'bold', size = 8, color = 'black'),
        axis.text.y = element_text(size = 7, colour = 'black'),
        axis.title.y = element_text(face = 'bold', size = 8),
        legend.title = element_text(face = 'bold', size = 8),
        legend.text = element_text(size = 8),
        legend.text.align = 0.5,
        legend.position = 'bottom',
        #legend.box.background = element_rect(colour = "black"),
        legend.title.align = 0.5,
        strip.text = element_text(face = 'bold', size = 8))
DivCor_plot

ggsave(plot = DivCor_plot, filename = 'Outputs/Figures/Fig3_DiversityCorrelations.pdf', width = 8.5, height = 8.5, units = 'cm', dpi=600)

# Plot boxplots showing variability within and between field dimensions
# Field_dim_plot <- ggplot(DivCor2, aes(x = factor(Field_dimension, levels = c('Tax', 'Func'), labels = c('Taxonomic', 'Functional')),
#                                       y = r)) +
#   geom_boxplot() +
#   #geom_point(size = 2) +
#   geom_jitter(width = 0.1) +
#   facet_wrap(~factor(Spec_region, level = c('VNIR','SWIR'))) +
#   scale_x_discrete("", labels = c("Taxonomic", 'Functional')) +
#   scale_y_continuous("Correlation with Spectral Diversity \n(Pearson's r)", limits = c(-0.25,0.75), breaks = c(-.25,0,.25, .5,.75,1))+
#   theme_classic() +
#   theme(axis.ticks.x=element_blank(),
#         axis.text.x=element_text(size=12, color='black', face = 'bold'),
#         axis.text.y = element_text(size=10, colour = 'black'),
#         axis.title=element_text(face='bold', size=12),
#         legend.title=element_text(face='bold', size = 12),
#         legend.text.align = 0.5,
#         legend.position = 'bottom',
#         #legend.box.background = element_rect(colour = "black"),
#         legend.title.align=0.5,
#         strip.text = element_text(face = 'bold', size = 12))
#Field_dim_plot
# ggsave(plot = Field_dim_plot, filename =  "Outputs/Figures/DiversityCorrelations_field.jpeg", width=4, height=4, dpi=600)

# Plot boxplots showing variability within and between spectral metrics
# Spec_met_plot <- ggplot(DivCor2, 
#                  aes(x = factor(Spec_metric, levels = c('CV', 'CHV','SV')), 
#                      y = r, 
#                      fill = factor(Spec_metric, levels = c('CV', 'CHV','SV'), labels = c('Less heavily weights extremes', 'More heavily weights extremes', '')))) +
#   geom_boxplot() +
#   geom_jitter(width = 0.1) +
#   facet_wrap(~factor(Spec_region, level = c('VNIR','SWIR'))) +
#   scale_x_discrete("", labels = c("CV", 'CHV', 'SV')) +
#   scale_y_continuous("Correlation with Field-based Diversity \n(Pearson's r)", limits = c(-0.25,0.75), breaks = c(-.25,0,.25, .5,.75,1))+
#   theme_classic() +
#   labs(fill = 'Metric Sensitivity to Extremes') +
#   scale_fill_manual(values = c('grey60','grey40', 'grey40')) +
#   guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
#   theme( 
#         axis.text.x=element_blank(),
#         axis.text.y = element_text(size=10, colour = 'black'),
#         axis.title=element_text(face='bold', size=12),
#         legend.title=element_text(face='bold', size = 12),
#         legend.text = element_text(size = 12),
#         legend.text.align = 0.5,
#         legend.position = 'bottom',
#         #legend.box.background = element_rect(colour = "black"),
#         legend.title.align=0.5,
#         strip.text = element_text(face = 'bold', size = 12))
# Spec_met_plot
# ggsave(plot = Spec_met_plot, filename = "Outputs/Figures/DiversitCorrelations_SpecMet.jpeg", width=4, height=4, dpi=600)

rm(DivCor, DivCor.Test, DivCor.Test2, DivCor.total, DivCor.Total2, DivCor2, DivCor3, DivCor2Sum, DivCorSum)

#### 1.2.2.1__Plot Relationships w VNIR-Diversity --------------------------------------------------------------------
#   Further explore relationships between VNIR-spectral diversity and taxonomic/functional diversity; Figure S6

# Centre and standardize all spectral diversity indices and plot spectral ~ field 
DivInd.stand<-decostand(DivIndices[,11:16], method = "standardize", na.rm = T) %>%
  cbind(plot_field_id) %>%
  left_join(SiteList) %>%
  left_join(DivIndices[,1:10])

DivInd.stand2 <- DivInd.stand %>%
  gather(Spec_metric, value_stand, CV_VNIR:SV_SWIR) %>%
  separate(Spec_metric, c("Spec_metric",'Spec_region'))

# This function creates Spectral Diversity ~ Field Diversity plots
DivPlots <- function(df, # Function written specifically for DivInd.stand2, won't work with differently structured dataframes.
                     metric, # Field diversity metric
                     spectral # Spectral region
                     ){
  
  Plot <- ggplot(data = df %>% dplyr::filter(Spec_region == {{spectral}}),
                 aes(y = value_stand, x = {{metric}})) +
    geom_point(aes(fill = Spec_metric, shape = Spec_metric)) +
    geom_smooth(aes(colour = Spec_metric, fill = Spec_metric),method = 'lm', fullrange = T, se = F) +
    scale_shape_manual(values = c(21,22,24),  name ='Spectral Metric') +
    scale_colour_manual(values = c("black","#636363","#bdbdbd"), name = 'Spectral Metric') +
    scale_fill_manual(values = c("black","#636363","#bdbdbd"), name = 'Spectral Metric') +
    scale_y_continuous(expand = c(0,0), breaks = seq(-2, 4, by = 1), limits = c(-2.5,4)) +
    ylab("") +
    guides(fill = guide_legend(title.position="top", title.hjust = 0.5), 
           colour = guide_legend(title.position="top", title.hjust = 0.5),
           shape =  guide_legend(title.position="top", title.hjust = 0.5)) +
    theme_classic() +
    theme(axis.title.x = element_markdown(face = 'bold', size = 12), 
          axis.text.x = element_text(size = 10, colour = 'black'),
          title = element_text(face = 'bold', size = 10), 
          axis.text.y = element_text(size = 10, colour = 'black'), 
          axis.title.y = element_text(face = 'bold', size = 12), 
          legend.title = element_text(face = 'bold', size = 10), 
          legend.position = 'none',
          #legend.box.background = element_rect(colour = "black"),
          legend.title.align = 0.5)
  
  return(Plot)
}

DivInd_VNIRPlots <- ggarrange(DivPlots(DivInd.stand2, q0, 'VNIR') +
                    labs(title = 'Species Richness', x = '<sup>0</sup>*D*') +
                    scale_x_continuous(expand = c(0,0), limits = c(0,15), breaks = seq(0,15, by =5)),
                  DivPlots(DivInd.stand2, q1, 'VNIR') +
                    labs(title = 'Exp. Shannon Index', x = '<sup>1</sup>*D*') +
                    scale_x_continuous(expand = c(0,0),limits = c(0,8), breaks = seq(0,8, by = 2)),
                  DivPlots(DivInd.stand2, q2, 'VNIR') +
                    labs(title = 'Inv. Simpson Index', x = '<sup>2</sup>*D*') +
                    scale_x_continuous(expand = c(0,0), limits = c(0,6), breaks = seq(0,6,by = 2)),
                  DivPlots(DivInd.stand2, J, 'VNIR') +
                    labs(title = 'Jaccards Evenness', x = "J\'") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = .25)),
                  DivPlots(DivInd.stand2, FRic, 'VNIR') +
                    labs(title = 'Functional Richness', x = "F<sub>Ric</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,33), breaks= seq(0,30, by = 10)),
                  DivPlots(DivInd.stand2, FEve, 'VNIR') +
                    labs(title = 'Functional Evenness', x = "F<sub>Eve</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by =0.25)),
                  DivPlots(DivInd.stand2, FDiv, 'VNIR') +
                    labs(title = 'Functional Divergence', x = "F<sub>Div</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0.25,1), breaks = seq(0.25,1, by = 0.25)),
                  DivPlots(DivInd.stand2, FDis, 'VNIR') +
                    labs(title = "Functional Dispersion", x ="F<sub>Dis</sub>") +
                    scale_x_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                  DivPlots(DivInd.stand2, RaoQ, 'VNIR') +
                    labs(title = 'Rao Entropy') +
                    scale_x_continuous(expand = c(0,0), limits = c(0,12), breaks = seq(0,12, by = 4)),
                  common.legend = TRUE, legend = 'bottom', labels = 'AUTO'
) 

DivInd_VNIRPlots <- annotate_figure(DivInd_VNIRPlots, left = text_grob("Standardized VNIR-Spectral Diversity", color = 'black', face = 'bold', size = 12, rot = 90))
DivInd_VNIRPlots

ggsave(plot = DivInd_VNIRPlots, filename = 'Outputs/Figures/FigS6_DiversityCorrelations_VNIR.jpeg', width= 7.5, height = 7.5, dpi=600)

rm(DivInd.stand, DivInd.stand2, DivIndPlots)

# 2.0__Environmental Drivers of Community Properties ---------------------------------------------
#   Compare whether biodiversity dimensions are driven by the same environmental gradients - AT MONT MEGANTIC ONLY
#       - RDAs for composition and GAMs for diversity.

# But first, examine the correlation between environmental variables
EnvrCor <- cor(Envr[,4:9], use="complete.obs") # Cols 4:9 contain environmental variables

corrplot(EnvrCor, method = 'number') # Roughness and slope are strongly correlated ... don't model roughness!

# Make sure environmental variables are dimensionally homogeneous (e.g., centered and standardized), so regression coefficients can be directly compared.
Envr.MMG <- Envr %>%
  dplyr::filter(site =='MtMeg-1')

Envr.MMG.stand <- decostand(Envr.MMG[,4:9], method = 'standardize') # Cols 4:9 contain environmental variables

## 2.1__Composition  --------------------------------------------------------
# Redundancy analysis (RDA) is a direct gradient analysis method that summarises linear relationships between components of response variables that are explained by a set of explanatory variables.
# The total variance of the data set is partitioned into constrained and unconstrained variances,
#   -If constrained variance > unconstrained variance then analysis suggests that variation in response variables is accounted by explanatory variables.
#   -If unconstrained variance > constrained variance then the results should be interpreted with caution.
# Significance value for the overall RDA solution and individual RDA axes can be determined by permutation.

# This function fits RDAs and creates dataframes with different summary statistics related to:  
#   a) the overall model,
#   b) the canonical axes, 
#   c) the explanatory variables


### 2.1.1__Fit RDAs ---------------------------------------------------------

# This function fits the Comp ~ environment model and extracts necessary model outputs, it is not coded to be generlizably (aka specific to envr variables measured here) 
RDA <- function(comp, # Composition dataframe 
                envr, # Envr dataframe
                site, # Character value stating what site the data is from (i.e., MMG)
                dimension # Biodiversity dimension of composition dataframe (i.e., tax, fun, or spectral)
                ){ 
  
  # a) run RDA model
  rda <- rda(comp ~ elevation + slope  + northness + eastness + TWI  , data = envr)
  
   # b) test model significance via permutations 
  rda.anova <- anova.cca(rda, permutations = 999)
  rda.anova.axes <- anova.cca(rda, by = 'axis', permutations = 999)
  rda.anova.mar <- anova(rda, by = 'mar', permuations = 999) 
  
  # c) create objects with model outputs
  rda.summary <- summary(rda)
  adj.R2 <- RsquareAdj(rda)$adj.r.squared
  fvalue <- rda.anova[1,3]
  pvalue <- rda.anova[1,4]
  constrained <- rda.summary[['constr.chi']] / rda.summary[['tot.chi']]
  unconstrained <- rda.summary[["unconst.chi"]] / rda.summary[['tot.chi']]
  
  # 1d) create dataframe with overall model statistics
  rda.output <- data.frame(site = site,
                           dimension = dimension,
                           adj.R2 = adj.R2,
                           fvalue = fvalue,
                           pvalue = pvalue,
                           constrained = constrained,
                           unconstrained = unconstrained)
  
  # 2d) create dataframe with axis model statistics
  rda.axes <- rda.anova.axes %>%
    rownames_to_column(var = 'axis') %>%
    mutate(site = {{site}}, .before = axis) %>%
    mutate(dimension = {{dimension}}, .before = axis)
  
  # 3d) create dataframe with explan. variable model statistics
  rda.explanvar <- rda.anova.mar %>% 
    rownames_to_column(var = 'explanvar') %>%
    mutate(site = {{site}}, .before = explanvar) %>%
    mutate(dimension = {{dimension}}, .before = explanvar)
  
  # 4d) create dataframe with how explan. variables are correlated with axes
  rda.correlations <- data.frame(rda[["CCA"]][['biplot']]) %>% 
    rownames_to_column(var = 'explanvar') %>%
    mutate(site = {{site}}, .before = RDA1) %>%
    mutate(dimension = {{dimension}}, .before = RDA1)
  
  
  return(list(rda.output, rda.axes, rda.explanvar,rda.correlations, 
              rda))
}


#### 2.1.1.1__Taxonomic Composition RDA ---------------------------------------------

# Only focused on MMG, too large of envr differences to combine MSB and MMG and MSB has too small sample size to be modelled alone
TaxComp.MMG <- TaxComp %>%
  dplyr::filter(site == 'MtMeg-1')

# Hellinger transform abundances and remove speces not observed at MMG
TaxComp.MMG.hellinger <-sqrt(TaxComp.MMG[,4:33]) %>% # Cols 4:33 contain species abundances
    dplyr::select(-c("Carya.cordiformis..Wangenheim..K..Koch", # This removes species that were not observed at Mont Megantic but were at Mont St Bruno
                     "Larix.laricina..Du.Roi..K..Koch",
                     "Pinus.strobus.Linnaeus",
                     "Quercus.rubra.Linnaeus",
                     "Tilia.americana.Linnaeus",
                     "Ulmus.rubra.Muhlenberg")) 

# Run RDA 
Tax.MMG.rda <- RDA(TaxComp.MMG.hellinger, Envr.MMG.stand, 'MMG', 'Taxonomic')

#### 2.1.1.2__Functional Composition RDA --------------------------------------------
# Join informatin about sites and filter for MMG
FunComp.MMG <- FunComp %>%
  left_join(SiteList) %>%
  dplyr::filter(site =='MtMeg-1')

# Standardize trait values
FunComp.MMG.stand <- decostand(FunComp.MMG[,2:16], method = 'standardize') 

# Run RDA
Fun.MMG.rda <- RDA(FunComp.MMG.stand, Envr.MMG.stand, 'MMG', 'Functional')

# Join output tables across dimensions
RDA.overall <- bind_rows(Tax.MMG.rda[[1]], Fun.MMG.rda[[1]])
RDA.axes <- bind_rows(Tax.MMG.rda[[2]], Fun.MMG.rda[[2]])
RDA.explanvar <- bind_rows(Tax.MMG.rda[[3]], Fun.MMG.rda[[3]])
RDA.correlations <- bind_rows(Tax.MMG.rda[[4]], Fun.MMG.rda[[4]])


#### 2.1.1.3__VNIR-Spectral Composition RDA --------------------------------------------------
# Join information about sites and filter for MMG
VNIRComp.MMG <- VNIRComp %>%
  left_join(SiteList) %>%
  dplyr::filter(site == 'MtMeg-1') %>%
  relocate(site, .after = plot_field_id)

# Run RDA
VNIR.MMG.rda <- RDA(VNIRComp.MMG[,-c(1:2)], # Don't include the cols that contain plot and site information 
                    Envr.MMG.stand, 'MMG', 'VNIR')

# Join output tables across dimensions
RDA.overall   <- bind_rows(RDA.overall, VNIR.MMG.rda[[1]])
RDA.axes      <- bind_rows(RDA.axes, VNIR.MMG.rda[[2]]) 
RDA.explanvar <- bind_rows(RDA.explanvar, VNIR.MMG.rda[[3]])
RDA.correlations <- bind_rows(RDA.correlations, VNIR.MMG.rda[[4]])

#### 2.1.1.4__SWIR-Spectral Composition RDA --------------------------------------------------
# Join information about sites and filter for MMG
SWIRComp.MMG <- SWIRComp %>%
  left_join(SiteList) %>%
  dplyr::filter(site == 'MtMeg-1') %>%
  relocate(site, .after = plot_field_id)

# Run RDA
SWIR.MMG.rda <- RDA(SWIRComp.MMG[,-c(1:2)], # Don't include the cols that contain plot and site information 
                    Envr.MMG.stand, 'MMG', 'SWIR')

# Join output tables across dimensions
RDA.overall <- bind_rows(RDA.overall, SWIR.MMG.rda[[1]])
RDA.axes      <- bind_rows(RDA.axes, SWIR.MMG.rda[[2]]) 
RDA.explanvar <- bind_rows(RDA.explanvar, SWIR.MMG.rda[[3]])
RDA.correlations <- bind_rows(RDA.correlations, SWIR.MMG.rda[[4]])

write.csv(RDA.overall, file = "Outputs/Statistics/TableS3_RDA_Stats_Overall.csv")
write.csv(RDA.axes, file = "Outputs/Statistics/RDA_Stats_Axes.csv")
write.csv(RDA.explanvar, file = 'Outputs/Statistics/TableS4_RDA_Stats_ExplanatoryVars.csv')
write.csv(RDA.correlations, file = 'Outputs/Statistics/RDA_Stats_Correlations.csv')

# Create wide format for Table 3
RDA_explanvar2 <- RDA_explanvar %>%
  select(c(dimension, explanvar, Variance)) %>%
  filter(explanvar != 'Residual') %>%
  group_by(dimension) %>%
  mutate(proportion = Variance/ sum(Variance) * 100) %>%
  mutate(proportion = round(proportion, 2)) %>%
  select(-Variance) %>%
  spread(key = "dimension", value = "proportion") %>%
  arrange(desc(Taxonomic)) %>%
  relocate(Taxonomic, .after = explanvar) %>%
  relocate(Functional, .after = Taxonomic)

write.csv(RDA_explanvar2, file = 'Outputs/Statistics/Table3_RDA_Stats_ExplanatoryVars.csv')


### 2.1.2__Plot RDAs ------------------------------------------------------
# Make RDA triplots for Taxonomic, Functional, and VNIR for MMG
# Determine the percentage of conifer abundance in each plot - to aide the interpretation of trioplots
PercConifer.MMG <- PercConifer %>%
  filter(site == "MtMeg-1")

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
windows(width = 3.34646, height = 3 * 3.34646)

RDAtriplots <- ggarrange(plotRDA2(plotRDA1(TaxComp.MMG.hellinger, Envr.MMG.stand)) +
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
                        
                        plotRDA2(plotRDA1(VNIRComp.MMG[,-c(1:2)], Envr.MMG.stand)) + 
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

ggsave(plot = RDAtriplots, filename = 'Outputs/Figures/Fig4_RDAtriplots_nolabels.pdf', width=3.346, height=9.44882, dpi=600)

rm(VNIRComp.MMG, VNIR.MMG.rda, EnvrCor, Fun.MMG.rda, FunComp.MMG, FunComp.MMG.stand, RDA.axes, 
   RDA.correlations, RDA.explanvar, RDA.overall, RDAtriplots, SWIRComp.MMG, SWIR.MMG.rda, TaxComp.MMG.hellinger, Tax.MMG.rda)


## 2.2__Diversity -----------------------------------------------------------
# Generalized Additive Models (GAMs) is a GLM in which the linear response variable depends linearly on unknown smooth functions of some predictor variables. 
# The relationship between individual predictors and the response variable follow smooth patterns that can be linear or non-linear. There are two pros:
#     i)  Flexible predictor functions can uncover hidden patterns in data
#     ii) Regularization of predictor functions helps avoid overfitting. 
# There are two ways of estimating the smoothing parameter:
#     i)  Generalized cross-validation (GCV)  
#     ii) Restricted maximum likelihood (REML) - Here, I will use REML because GCV is prone to under-smoothing
# Note: package mgcv uses penalized regression splines and by default defines basis functions 


### 2.2.1__Fit GAMs ---------------------------------------------------------
#DivIndices <- read.csv("Outputs/Data/DiversityIndices.csv")

# Subset data for just Mont Megantic and merge the diversity dataframe with the standardized environmental variable dataframe
DivIndices.MMG <- DivIndices %>%
  left_join(SiteList) %>%
  dplyr::filter(site == 'MtMeg-1') %>%
  cbind(Envr.MMG.stand) %>%
  left_join(PercConifer.MMG)


##### 2.2.1.1__Taxonomic Diversity GAM ----------------------------------------
# Taxonomic Diversity (q1)
# Fit model
GAM.q1  <- gam(q1 ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.q1)

# Check model assumptions
Check.GAM.q1 <- getViz(GAM.q1) 
check(Check.GAM.q1,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.q1)

# Visualize smoothed fit
plot(GAM.q1, residuals = T, pch = 19, cex=.3)


#### 2.2.1.2__Functional Diversity GAM ---------------------------------------
# Functional Diversity (FDis)
# Fit model
GAM.FDis  <- gam(FDis ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.FDis)

# Check assumptions
Check.GAM.FDis <- getViz(GAM.FDis) 
check(Check.GAM.FDis,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.FDis)

# Visualize smoothed fit
plot.gam(GAM.FDis, residuals = T, pch = 19, cex=.3)

##### 2.2.1.3__VNIR-Spectral Diversity GAM -------------------------------------
# VNIR Diversity (SV_VNIR)
# Fit model
GAM.SV.VNIR  <- gam(SV_VNIR ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.SV.VNIR)

# Check assumptions
Check.GAM.SV.VNIR <- getViz(GAM.SV.VNIR) 
check(Check.GAM.SV.VNIR,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.SV.VNIR)

# Visualize smoothed fit
plot.gam(GAM.SV.VNIR, residuals = T, pch = 19, cex=.3)

# Visualize smoothed model
plot(GAM.SV.VNIR, residuals = T, pch = 19, cex=.3)

#### 2.2.1.4__SWIR-Spectral Diversity GAM ------------------------------------
# SWIR Diversity (SV_SWIR)
# Fit model
GAM.SV.SWIR  <- gam(SV_SWIR ~ s(elevation) + TWI + northness + eastness + slope, data = DivIndices.MMG, method = 'REML')
summary(GAM.SV.SWIR)

# Check assumptions
Check.GAM.SV.SWIR <- mgcViz::getViz(GAM.SV.SWIR) 
check(Check.GAM.SV.SWIR,
      a.qq = list(method = "tnorm", 
                  a.cipoly = list(fill = "light blue")), 
      a.respoi = list(size = 0.5), 
      a.hist = list(bins = 10))

gam.check(GAM.SV.SWIR)

# Visualize smoothed fit
plot(GAM.SV.SWIR, residuals = T, pch = 19, cex=.3)

# This function fits GAMs and extracts GAM model outputs
GAMextract <- function(df, # Dataframe containing diversity and environmental variables per plot
                       metric, # The diversity metric you want to model
                       site, # A character value stating what site the data is from
                       dimension # A character value stating the dimension of diversity index (i.e., tax, fun, spectral)
                       ){
  df2 <- df %>%
    rename(index = metric)
  
  # 1a) Fit full model 
  full <- gam(index ~ s(elevation) + TWI + northness + eastness + slope, method = 'REML', data = df2)
  
  # 2a) Fit reduced models, but keep smoothing parameter (sp) from the full model 
  elevation <- gam(index ~ TWI + northness + eastness + slope, method = 'REML', data = df2)
  TWI <- gam(index ~ s(elevation) + northness + eastness + slope, sp = c(full$sp[1]), method = 'REML', data = df2)
  northness <- gam(index ~ s(elevation) + TWI + eastness + slope, sp = c(full$sp[1]), method = 'REML', data = df2)
  eastness <- gam(index ~ s(elevation) + TWI + northness + slope, sp = c(full$sp[1]), method = 'REML', data = df2)
  slope <- gam(index ~ s(elevation) + TWI + northness + eastness, sp = c(full$sp[1]), method = 'REML', data = df2)
  null <- gam(index ~ 1, data = df2)

  # b) Calculate the proportion of deviance explained by each variable
  devElev <- (elevation$deviance - full$deviance) / null$deviance
  devTWI  <- (TWI$deviance - full$deviance) / null$deviance
  devNorthness <- (northness$deviance - full$deviance) / null$deviance
  devEastness <- (eastness$deviance - full$deviance) / null$deviance
  devSlope <- (slope$deviance - full$deviance) / null$deviance
  
  # 1c) Extract model outputs
  modelsummary <- summary(full)
  adj.R2 <- modelsummary$r.sq
  dev.expl <- modelsummary$dev.expl
  AIC <- AIC(full)
  GCV <- full[["gcv.ubre"]][["REML"]]
  
  # 2c) Create dataframe with overall model statistics
  GAMoutput <- data.frame(site = site,
                           dimension = dimension,
                           index = metric,
                           adj.R2 = adj.R2,
                           dev.expl = dev.expl,
                           AIC = AIC,
                           GCV = GCV)
  
  # d) Create dataframe with explan. variable statistics
  GAMexplanvar <- data.frame(metric = metric,
                             explanvar = c('Elevation', 'TWI', 'Northness', 'Eastness', 'Slope'),
                             deviance = c(devElev, devTWI, devNorthness, devEastness, devSlope)) %>%
    mutate(proportion = deviance/sum(deviance) * 100) # Proportional deviance explained to compare across models

  return(list(GAMoutput, GAMexplanvar))
}

# Fit GAMs for q1, FDis, and SV
Tax.MMG.gam <- GAMextract(DivIndices.MMG, 'q1', 'MMG', 'Taxonomic')
Fun.MMG.gam <- GAMextract(DivIndices.MMG, 'FDis', 'MMG', 'Functional')
VNIR.MMG.gam <- GAMextract(DivIndices.MMG, 'SV_VNIR', 'MMG', 'VNIR')
SWIR.MMG.gam <- GAMextract(DivIndices.MMG, 'SV_SWIR', 'MMG', 'SWIR')

GAM_overall <- bind_rows(Tax.MMG.gam[[1]], Fun.MMG.gam[[1]], VNIR.MMG.gam[[1]], SWIR.MMG.gam[[1]])

GAM_explanvar <- bind_rows(Tax.MMG.gam[[2]], Fun.MMG.gam[[2]], VNIR.MMG.gam[[2]], SWIR.MMG.gam[[2]])

# Create wide format for Table 4
GAM_explanvar2 <- GAM_explanvar %>%
  select(c(metric, explanvar, proportion)) %>%
  mutate(proportion = round(proportion, 2)) %>%
  spread(key = "metric", value = "proportion") %>%
  arrange(desc(q1)) %>%
  relocate(q1, .before = FDis) %>%
  relocate(SV_VNIR, .after = FDis)
  
write.csv(GAM_overall, file = "Outputs/Statistics/GAMS_Overall.csv", row.names = F)
write.csv(GAM_explanvar, file = "Outputs/Statistics/TableS5_GAM_Explanvar.csv", row.names = F)
write.csv(GAM_explanvar2, file = 'Outputs/Statistics/Table4_GAM_Explanvar.csv', row.names = F)

### 2.2.2__Plot GAMs ---------------------------------------------------------

# Create a dataframe with unstandardized (raw) environmental variables
DivIndices.MMG.raw <- DivIndices %>%
  left_join(SiteList) %>%
  dplyr::filter(site == 'MtMeg-1') %>%
  merge(Envr.MMG) %>%
  left_join(PercConifer.MMG) 

# This function plots the relationship with elevation (i.e., the smooth variable in the GAMs)
plotGAM <- function(metric, # Diversity metric of interest (x-axis)
                    envrvar, # Environmental variable of interest (y-axis)
                    df, # Dataframe containing diversity and standardized environmental variables per plot - for fitting models
                    df3, # Dataframe containing diversity and raw environmental variables per plot - for plotting
                    smoothed # If smoothed is TRUE, it will plot the smoothed function fit by GAM, if FALSE then it will plot the linear relationship.
                    ){
  
  # Rename metric, to be able to fit the model
  df2 <- df %>%
    rename(index = {{metric}})
  
  # Fit the diversity-environment model
  full <- gam(index ~ s(elevation) + TWI + northness + eastness + slope, method = 'REML', data = df2)
  
  # This creates a new object that is a character value of envr variable of interest
  envrvar2 <- deparse(substitute(envrvar))
  
  # Create object of the smoothing parameters of each of the envr predictors
  smooth <- ifelse(envrvar2 == 'elevation', full$sp[1],
                   ifelse(envrvar2 == 'TWI', full$sp[2],
                          ifelse(envrvar2 == 'northness', full$sp[3],
                                 ifelse(envrvar2 == 'eastness', full$sp[4], full$sp[5]))))
  
  if(smoothed == TRUE){ 
    
  # Plot diversity-envr fit with raw data behind
  Plot <- ggplot(df3, aes(x = {{envrvar}}, y = {{metric}}))+
    stat_smooth(method = gam, formula=y~s(x, sp = smooth), se=T,colour='black') +
    #stat_smooth(method = lm, se=T, colour='black') +
    geom_point(aes(fill = PercCon), size = 2, pch = 21) +
    scale_fill_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
    #scale_y_continuous(expand=c(0,0)) + 
    #scale_x_continuous(expand=c(0,0)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(size = 8, color = 'black'), 
          axis.text.y = element_text(size = 8, color = 'black'),
          axis.title.y = element_markdown(face = 'bold', size = 10),
          axis.title.x = element_text(face = 'bold', size = 10),
          legend.title = element_text(face = 'bold', size = 10),
          legend.text = element_text(face = 'bold', size= 10, color ='black'),
          plot.title = element_text(face = "bold"),
          legend.position = 'none') 
  }
  
  if(smoothed == FALSE){
    
  Plot <- ggplot(df3, aes(x={{envrvar}}, y={{metric}}))+
    #stat_smooth(method=gam, formula=y~s(x), se=T,colour='black') +
    stat_smooth(method = lm, se=T, colour='black') +
    geom_point(aes(fill = PercCon), size = 3, pch = 21) +
    scale_fill_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, barwidth = 10)) +
    #scale_y_continuous(expand=c(0,0)) + 
    #scale_x_continuous(expand=c(0,0)) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(size=8, color='black'),
          axis.text.y = element_text(size=8, color='black'),
          axis.title.y = element_markdown(face='bold', size=10),
          axis.title.x = element_text(face='bold', size=10),
          legend.title = element_text(face='bold', size = 10),
          legend.text = element_text(face = 'bold', size= 10, color ='black'),
          plot.title = element_text(face = "bold"),
          legend.position = 'none') 
  }
  return(Plot)
}

# Plot diversity-elevation relationships, Figure 5.
GAMElevation <- ggarrange(plotGAM(q1, elevation, DivIndices.MMG, DivIndices.MMG.raw, TRUE) +
                            labs(title = 'Taxonomic', x='Elevation (m a.s.l.)', y='Exp. Shannon Index (<sup>1</sup>*D*)') +
                            scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                            scale_y_continuous(expand = c(0,0), limits = c(0,7), breaks = seq(0,7, by = 1)), 
                          
                          plotGAM(FDis, elevation, DivIndices.MMG, DivIndices.MMG.raw, TRUE)+
                            labs(title = 'Functional',  x='Elevation (m a.s.l.)', y='Functional Dispersion (F<sub>Dis</sub>)') +
                            scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                            scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                          
                          plotGAM(SV_VNIR, elevation, DivIndices.MMG, DivIndices.MMG.raw, TRUE)+
                            labs(title = 'VNIR-Spectral', x='Elevation (m a.s.l.)', y='Spectral Variance (SV)') +
                            scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) ,
                          #scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                          
                          nrow = 3, ncol = 1,
                          common.legend = T, legend = 'bottom',
                          labels =  'AUTO')

GAMElevation

ggsave(plot = GAMElevation, filename = 'Outputs/Figures/Fig5_GAM_plots_elevation.pdf', width=3.346, height=9.44882, dpi=600)

# Plot diversity-environment relationships, Figure S7.
GAMPlot <- ggarrange(plotGAM(q1, elevation, DivIndices.MMG, DivIndices.MMG.raw, TRUE) +
                       labs(title = 'A. Taxonomic', x='', y='Exp. Shannon Index (<sup>1</sup>*D*)') +
                       scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)) ,
                     plotGAM(q1, northness, DivIndices.MMG, DivIndices.MMG.raw, FALSE) +
                       labs(title = '', x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)) ,
                     plotGAM(q1, eastness, DivIndices.MMG, DivIndices.MMG.raw, FALSE) +
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)),
                     plotGAM(q1, slope, DivIndices.MMG, DivIndices.MMG.raw, FALSE) +
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(0,30), breaks = seq(0,30, by = 5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)),
                     plotGAM(q1, TWI, DivIndices.MMG, DivIndices.MMG.raw, FALSE) +
                       labs(title = '', x='', y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(2,12), breaks = seq(0,12, by = 2)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8, by = 1)),

                     plotGAM(FDis, elevation,  DivIndices.MMG, DivIndices.MMG.raw, TRUE)+
                       labs(title = 'B. Functional', x='', y='Functional Dispersion (F<sub>Dis</sub>)') +
                        scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                        scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM(FDis, northness, DivIndices.MMG, DivIndices.MMG.raw, FALSE)+
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM(FDis, eastness,  DivIndices.MMG, DivIndices.MMG.raw, FALSE)+
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM(FDis, slope, DivIndices.MMG, DivIndices.MMG.raw, FALSE) +
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(0,30), breaks = seq(0,30, by = 5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),
                     plotGAM(FDis, TWI, DivIndices.MMG, DivIndices.MMG.raw, FALSE)+
                       labs(title = '',x='',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(2,12), breaks = seq(0,12, by = 2)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4, by = 1)),

                     plotGAM(SV_VNIR, elevation,  DivIndices.MMG, DivIndices.MMG.raw, TRUE)+
                       labs(title = 'C. VNIR-Spectral', x='Elevation\n(m a.s.l.)', y='Spectral Variance (SV)') +
                       scale_x_continuous(expand = c(0,0), limits = c(400,1100), breaks = seq(400, 1000, by = 200)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM(SV_VNIR, northness,  DivIndices.MMG, DivIndices.MMG.raw, FALSE)+
                       labs(title = '',x='Northness\n(cos(Aspect))',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM(SV_VNIR, eastness,  DivIndices.MMG, DivIndices.MMG.raw, FALSE)+
                       labs(title = '',x='Eastness\n(sin(Aspect))',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(-1,1), breaks = seq(-1,1, by = 0.5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM(SV_VNIR, slope,  DivIndices.MMG, DivIndices.MMG.raw, FALSE) +
                       labs(title = '',x='Slope\n(\u00B0)',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(0,30), breaks = seq(0,30, by = 5)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),
                     plotGAM(SV_VNIR, TWI,  DivIndices.MMG, DivIndices.MMG.raw, FALSE)+
                       labs(title = '',x='Topographic Wetness Index\n',y='') +
                       scale_x_continuous(expand = c(0,0), limits = c(2,12), breaks = seq(0,12, by = 2)) +
                       scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1, by = 0.25)),

                     nrow = 3, ncol = 5,
                     common.legend = T, legend = 'bottom'
                     )

GAMPlot

ggsave(plot = GAMPlot, filename = 'Outputs/Figures/FigS7_GAM_plots.jpeg', width=14, height=11, dpi=300)


rm(Check.GAM.q1, Check.GAM.FDis, Check.GAM.SV.VNIR, Check.GAM.SV.SWIR, Tax.MMG.gam, Fun.MMG.gam, VNIR.MMG.gam, SWIR.MMG.gam, 
   GAM_explanvar, GAM_overall, GAM.FDis, GAM_explanvar2, GAM.q1, GAM.SV.VNIR, GAM.SV.SWIR, GAMElevation, GAMPlot)

# 3.0__Additional Figures and Tables ------------------------------------------------------

## 3.1__Study Site Figure -------------------------------------------------------

# Plot the Spectral Composition for Plot 1, Fig. 1d
Plot1_VNIR <- Spec_VNIR %>%
  dplyr::filter(plot_field_id == 'Plot01') %>%
  dplyr::select(c('454.2':'1059.08')) %>%
  gather(Wavelength, BandDepth) %>%
  group_by(Wavelength) %>%
  summarise(BDavg = mean(BandDepth),
            #BDsd = sd(BandDepth),
            BDmax = max(BandDepth),
            BDmin = min(BandDepth)) %>%
  mutate(Sensor = 'Casi',
         Wavelength = as.numeric(Wavelength))

Plot1_SWIR <- Spec_SWIR %>%
  dplyr::filter(plot_field_id == 'Plot01') %>%
  dplyr::select(c('972.5':'2412.5')) %>%
  gather(Wavelength, BandDepth) %>%
  dplyr::filter(Wavelength > 1059) %>%
  group_by(Wavelength) %>%
  summarise(BDavg = mean(BandDepth),
            #BDsd = sd(BandDepth),
            BDmax = max(BandDepth),
            BDmin = min(BandDepth)) %>%
  mutate(Sensor = 'Sasi',
          Wavelength = as.numeric(Wavelength))

SpectraPlot1 <- rbind(Plot1_VNIR, Plot1_SWIR)

SpectrumPlot <- ggplot(SpectraPlot1) +
                  geom_line(aes(x = Wavelength, y = BDavg), cex = 1.1) +
                  geom_ribbon(aes(x = Wavelength, ymin = BDmin, ymax = BDmax), alpha = 0.2) +
                  annotate("rect", xmin = 899, xmax = 957, ymin = 0, ymax = 1, colour = 'white', fill = 'white') + # These wavelengths are masked because noisy
                  annotate("rect", xmin = 1345, xmax = 1460, ymin = 0, ymax = 1, colour = 'white', fill = 'white') + # These wavelengths are masked because of water absorption
                  annotate("rect", xmin = 1790, xmax = 1950, ymin = 0, ymax = 1, colour = 'white', fill = 'white') + # These wavelengths are masked because of water absorptuib
                  annotate("segment", x = 454, xend = 1059, y = -0.07, yend = -0.07, cex = 0.5) +
                  annotate("text", x = 756.5, y = -0.06, vjust = 0, label = 'CASI', fontface =2, size = 5) +
                  annotate("segment", x = 972, xend = 2412, y = -0.08, yend = -0.08, cex = 0.5) +
                  annotate("text", x = 1692, y = -0.06, vjust = 0, label = 'SASI', fontface =2, size =5) +
                  scale_y_continuous(limits = c(-0.1,1), expand = c(0,0)) +
                  scale_x_continuous(limits = c(300,2500), expand = c(0,0), breaks = seq(400,2400,400)) +
                  theme_classic() +
                  xlab("Wavelength (nm)") +
                  ylab("Normalized Reflectance") +
                  theme(axis.title.x=element_text(face ='bold', size = 14),
                        axis.text.x = element_text(size = 12, colour = 'black'),
                        axis.text.y=element_text(size = 12, colour = 'black'),
                        axis.title.y =element_text(face = 'bold', size = 14))
SpectrumPlot

ggsave(plot = SpectrumPlot, filename = 'Outputs/Figures/Fig1d_Plot1_Spectrum.jpeg', width=5, height=4, dpi=600)

## 3.2__Northness, Elevation, and Red Spruce -------------------------------

# Plot northness-elevation relationship, Fig.S8a
Northness <- full_join(TaxComp.MMG, PercConifer.MMG) %>%
  merge(Envr.MMG) %>%
  mutate(elevation = round(elevation, 0))

Northness_plot <- ggplot(Northness, aes(x = northness, y = elevation, fill = PercCon, label = elevation))+
  geom_point(size = 2, pch = 21) +
  scale_fill_gradient2(low = 'white', mid = 'black', high = 'white', midpoint = 0.5, name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
  theme_classic() +
  ylab("Elevation (m a.s.l.) ")+
  xlab("Northness\n (cos(aspect))")+
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10)) +
  labs(title = "") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = 'black'),
        axis.text.y = element_text(size = 9, colour = 'black'),
        axis.title = element_text(face = 'bold', size = 10),
        legend.text = element_text(face ='bold', size = 10),
        legend.title = element_text(face ='bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5) 
Northness_plot

Northness_plot2 <- ggplot(Northness, aes(x = northness, y = Picea.rubens.Sargent, label = elevation, fill = elevation)) +
  geom_point(size = 2, pch = 21) +
  scale_fill_gradient2(low = 'white', mid = 'black', high = 'white', midpoint = median(Envr.MMG$elevation), name = "Elevation (m a.s.l)") +
  theme_classic() +
  ylab("Relative Abundance\n of Picea rubens (%) ")+
  xlab("Northness\n (cos(aspect))")+
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, barwidth = 10)) +
  labs(title = "") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = 'black'),
        axis.text.y = element_text(size = 8, colour = 'black'),
        axis.title = element_text(face = 'bold', size = 10),
        legend.text = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        legend.position = 'bottom',
        legend.title.align = 0.5) 
Northness_plot2

windows(width = 3.334646, height = 7)
Northness_plot3 <- ggarrange(Northness_plot,
                        Northness_plot2,
                        nrow = 2, ncol = 1,
                        common.legend = F, labels = "AUTO")
Northness_plot3

ggsave(plot = Northness_plot3, filename ='Outputs/Figures/FigS8_Northness.jpeg', width=3.334646, height=7, dpi=600)


## 3.3__Wavelength Contributions to Spec Div -------------------------
#   Visualize wavelengths contributing to SV, Figure S4

# This function calculates the variance in normalized reflectance at each band per plot
calc.wavVar <- function(df, # Dataframe containing normalized spectral reflectance
                        grouping_var, # Column within df that contains the variable that you want to group over, here it is plots
                        .x # Columns within df that contain wavelengths/bands at which reflectance was quantified
                        ){
  # 1) Calculate how many spectral points (pixels) are within the areas of interest (ie. plots)
  points <- df %>%
    group_by({{grouping_var}}) %>%
    summarise(points = n())
  
  # 2) Calculate the variance for each wavelength/band
  SV.df  <- df %>% group_by({{grouping_var}}) %>%
    summarise_all(~sum((.x - mean(.x))^2)) %>%
    rowwise({{grouping_var}}) %>%
    left_join(points) %>%
    summarise_all(~sum(.x / (points - 1)))
  
  return(SV.df)
}

# Calculate variance in normalized reflectance for each wavelength
#   a) VNIR wavelengths
SV_VNIR <- Spec_VNIR %>%
  dplyr::select(c(plot_field_id, '454.2':'1059.08')) %>%
  calc.wavVar(plot_field_id, c('454.2':'1059.08'))

SV_VNIR2 <- SV_VNIR %>%
  gather(Wavelength, Variance, c('454.2':'1059.08')) %>%
  group_by(Wavelength) %>%
  mutate(Varavg = mean(Variance)) %>%
  mutate(Wavelength = as.numeric(Wavelength)) %>%
  left_join(PercConifer) 

SV_plot_VNIR <- ggplot(SV_VNIR2) +
  geom_line(aes(x=Wavelength, y=Variance, group = plot_field_id, colour = PercCon)) +
  scale_colour_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
  geom_line(aes(x=Wavelength, y=Varavg), colour = 'black', linewidth = 1.1) +
  annotate("rect", xmin = 899, xmax = 957, ymin = -0.001, ymax = 0.15, fill= 'grey90') + # These wavelengths are masked
  ggbreak::scale_y_cut(breaks = 0.021, scales = 0.5, space  = .1) +
  scale_y_continuous(limits = c(-0.001, 0.15), expand = c(0, 0), breaks = c(0.0, 0.005, 0.01, 0.015, 0.02, .05, .10,.15)) +
  scale_x_continuous(limits = c(454, 1059), breaks = seq(500, 1000, 100)) +
  theme_classic() +
  xlab("Wavelength (nm)") +
  ylab("Variance") +
  labs(title = 'VNIR-Spectral') +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        plot.title = element_text(face = 'bold'),
        legend.position = 'none') 
SV_plot_VNIR

ggsave(plot = SV_plot_VNIR, filename = 'Outputs/Figures/FigS4a_Wavelenght_variance_VNIR.jpeg', width = 3.334646, height = 2.5, dpi=600)

#   b) SWIR wavelengths
SV_SWIR <- Spec_SWIR %>%
  select(c(plot_field_id, '972.5':'2412.5')) %>%
  calc.wavVar(plot_field_id, c('972.5':'2412.5'))

SV_SWIR2 <- SV_SWIR %>%
  gather(Wavelength, Variance, c('972.5':'2412.5')) %>%
  group_by(Wavelength) %>%
  mutate(Varavg = mean(Variance, na.rm = T)) %>%
  mutate(Wavelength = as.numeric(Wavelength)) %>%
  left_join(PercConifer)

SV_plot_SWIR <- ggplot(SV_SWIR2) +
  geom_line(aes(x=Wavelength, y=Variance, group = plot_field_id, colour = PercCon)) +
  scale_colour_gradient(low = '#bbe7a8', high = '#004a00', name = "Percent Coniferous Cover", limits = c(0,1), labels = c('0','','50','','100')) +
  geom_line(aes(x=Wavelength, y=Varavg), colour = 'black', linewidth = 1.1) +
  annotate("rect", xmin = 1345, xmax = 1460, ymin = -0.001, ymax = 0.15, fill= 'grey90') + # These wavelengths are masked
  annotate("rect", xmin = 1790, xmax = 1950, ymin = -0.001, ymax = 0.15, fill= 'grey90') + # These wavelengths are masked
  ggbreak::scale_y_cut(breaks = 0.021, scales = 0.5, space  = .1) +
  scale_y_continuous(limits = c(-0.001,0.15), expand = c(0,0), breaks =c(0.0, 0.005, 0.01, 0.015, 0.02, .05, .10,.15)) +
  scale_x_continuous(limits = c(1002, 2412), breaks = seq(1000, 2400, 400)) +
  theme_classic() +
  xlab("Wavelength (nm)") +
  ylab("Variance") +
  labs(title = 'SWIR-Spectral') +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        axis.title.x = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        plot.title = element_text(face = 'bold'),
        legend.position = 'none') 
SV_plot_SWIR

ggsave(plot = SV_plot_SWIR, filename = 'Outputs/Figures/FigS4b_Wavelenght_variance_SWIR.jpeg', width = 3.334646, height = 2.5, dpi = 600)


## 3.4__Species Mean Foliar Traits -----------------------------------------

# sp_means <- read.csv('~/PhD/Ch1/Trait_species_updated.csv')
# sp_sd <- read.csv('~/PhD/Ch1/Trait_species_updated_sd.csv')
# 
# traits_struct <- sp_means %>% 
#   mutate(SLA = paste(round(sp_means$trait_specific_leaf_area_m2_kg, 2), "+/-", round(sp_sd$trait_specific_leaf_area_m2_kg, 2))) %>%
#   mutate(LDMC = paste(round(sp_means$trait_leaf_dry_matter_content_mg_g, 2), "+/-", round(sp_sd$trait_leaf_dry_matter_content_mg_g, 2))) %>%
#   mutate(RWC = paste(round(sp_means$trait_leaf_relative_water_content_perc, 2), "+/-", round(sp_sd$trait_leaf_relative_water_content_perc, 2))) %>%
#   mutate(EWT = paste(round(sp_means$trait_equivalent_water_thickness_cm, 3), "+/-", round(sp_sd$trait_equivalent_water_thickness_cm, 3))) %>%
#   select(c('scientific_name', 'count_samples_leaf', 'SLA', 'LDMC', 'RWC','EWT'))
# write.csv(traits_struct, file = "Outputs/Data/TableS3_traits_struct.csv", row.names = F)
# 
# traits_C <- sp_means %>%
#   mutate(C = paste(round(sp_means$trait_c_perc, 2), "+/-", round(sp_sd$trait_c_perc, 2))) %>%
#   mutate(N = paste(round(sp_means$trait_n_perc, 2), "+/-", round(sp_sd$trait_n_perc, 2))) %>%
#   mutate(Sol = paste(round(sp_means$trait_soluble_perc, 2), "+/-", round(sp_sd$trait_soluble_perc, 2))) %>%
#   mutate(Hemi = paste(round(sp_means$trait_hemicellulose_perc, 2), "+/-", round(sp_sd$trait_hemicellulose_perc, 2))) %>%
#   mutate(Cel = paste(round(sp_means$trait_cellulose_perc, 2), "+/-", round(sp_sd$trait_cellulose_perc, 2))) %>%
#   mutate(Lig = paste(round(sp_means$trait_lignin_perc, 2), "+/-", round(sp_sd$trait_lignin_perc, 2))) %>%
#   mutate(Recal = paste(round(sp_means$trait_recalcitrants_perc, 2), "+/-", round(sp_sd$trait_recalcitrants_perc, 2))) %>%
#   select(c('scientific_name', 'count_samples_CN', 'C', 'N', 'Sol', 'Hemi','Cel', 'Lig', 'Recal'))
# write.csv(traits_C, file = "Outputs/Data/TableS4_traits_c.csv", row.names = F)
# 
# traits_pigs <- sp_means %>%
#   mutate(Chla = paste(round(sp_means$trait_chla_mg_g_disk_mass, 2), "+/-", round(sp_sd$trait_chla_mg_g_disk_mass, 2))) %>%
#   mutate(Chlb = paste(round(sp_means$trait_chlb_mg_g_disk_mass, 2), "+/-", round(sp_sd$trait_chlb_mg_g_disk_mass, 2))) %>%
#   mutate(Chlab = paste(round(sp_means$trait_chl_a_chl_b_ratio, 2), "+/-", round(sp_sd$trait_chl_a_chl_b_ratio, 2))) %>%
#   mutate(Carot = paste(round(sp_means$trait_carot_mg_g_disk_mass, 2), "+/-", round(sp_sd$trait_carot_mg_g_disk_mass, 2))) %>%
#   select(c('scientific_name', 'count_samples_pig', 'Chla', 'Chlb', 'Chlab', 'Carot'))
# write.csv(traits_pigs, file = "Outputs/Data/TableS5_traits_pigs.csv", row.names = F)


# 4.0_Full Range Spectra --------------------------------------------------
## 4.1__Full Spectral Composition ----------------------------------------

# Spectral composition is defined as the mean normalized reflectance across all bands (ie. wavelengths) per plot
VNIRComp2 <- VNIRComp %>%
  select(c(plot_field_id, '454.2':'970.77')) # Remove wavelengths that overlap with SASI sensor, 973.16 - 1059.08 nm

# Join VNIR and SWIR composition to create full range spectrum
FullComp <- left_join(VNIRComp2, SWIRComp) 

FullComp <- FullComp[match(plot_field_id, FullComp$plot_field_id), ]

# Change to long format and visualize full range spectrum
FullComp_long <- FullComp %>%
  dplyr::filter(plot_field_id == 'Plot01') %>%
  dplyr::select(-c('plot_field_id')) %>%
  gather(Wavelength, BandDepth) %>%
  mutate(Wavelength = as.numeric(Wavelength))

SpectrumPlot <- ggplot(FullComp_long) +
  geom_line(aes(x=Wavelength, y=BandDepth), cex=1.1) +
  annotate("rect", xmin = 899, xmax = 957, ymin = 0, ymax = 1, color="grey80", fill= 'grey80') +
  annotate("rect", xmin = 1345, xmax = 1460, ymin = 0, ymax = 1, color="grey80", fill= 'grey80') +
  annotate("rect", xmin = 1790, xmax = 1950, ymin = 0, ymax = 1, color="grey80", fill= 'grey80') +
  #annotate("rect", xmin = 972, xmax = 2435, ymin = -0.08, ymax = 1, alpha = .3, color="grey40", fill= 'grey20', linetype = 'dashed') +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
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

# Explore ordination
#           - First two PCs explain 97.67% of the variance
#           - PC1 69.92%, PC2 13.30%, PC3 10.50% 
#           - >99% variance explained by the PC1-PC8
pca.Full <- rda(FullComp[,-1])
summary(pca.Full) 

#   Extract wavelength scores and plot them along PCs
Full.band.scores <- data.frame(pca.Full$CA$v) %>%
  select(c(1:8)) %>%
  rownames_to_column(var="Bands") %>%
  mutate(Wavelength = gsub('BD_Band', "", Bands),
         Wavelength = gsub('nm', '', Wavelength),
         Wavelength = as.numeric(Wavelength))

PC1.Full <- ggplot(Full.band.scores, aes(x = Wavelength, y = PC1)) +
  geom_line() +
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC1.Full 

PC2.Full <- ggplot(Full.band.scores, aes(x = Wavelength, y = PC2)) +
  geom_line() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
PC2.Full 

### 4.1.1__Run Procrustes Analyses -----------------------------------------

# Load Tax, Fun, VNIR, and SWIR PCA outputs
PCA_outputs <- read.csv("Outputs/Data/PCA_outputs_updated.csv")

# Create dataframes of site scores across PCs 1-8
# PCs 1-8 contains greater than 90% of variation in composition across all dimensions 
FullPCA_scores <- PCA_sitescores(FullComp[,-1], plot_field_id, 'Full')

# Add full range PCA scores to other dimensions
PCA_outputs_full <- full_join(PCA_outputs, FullPCA_scores)

write.csv(PCA_outputs_full, "Outputs/Data/PCA_outputs_full.csv")


# Run procrustes analyses and procrustean permutation test (n = 999)
Tax.Full.out <- # Output of procrustes analysis btwn taxonomic and Full-spectral composition
  protest(
    X = PCA_outputs_full %>% # Tax site scores for PCs 1-8
      filter(dimension == 'Taxonomic') %>% 
      select(c(PC1:PC8))
    , 
    Y = FullPCA_scores[, 1:8], # VNIR-spectral site scores for PCs 1-8
    symmetric = T, 
    scores = "sites",
    permutations = 999
  )

Fun.Full.out <- # Output of procrustes analysis btwn functional and Full-spectral composition
  protest(
    X = PCA_outputs_full %>% # Fun site scores for PCs 1-8
      filter(dimension == 'Functional') %>% 
      select(c(PC1:PC8))
    , 
    Y = FullPCA_scores[, 1:8], # VNIR-spectral site scores for PCs 1-8
    symmetric = T, 
    scores = "sites",
    permutations = 999
  )

# Create a dataframe detailing procrustes output and merge with procustes outputs from other spectral comparisons
Procrustes_full <- data.frame(Field = c('Taxonomic', 'Functional'),
                         Spectral = c('Full', 'Full'),
                         r =c(Tax.Full.out$t0, Fun.Full.out$t0),
                         m12.squared =c(Tax.Full.out[["ss"]], Fun.Full.out[["ss"]]),
                         p = c(Tax.Full.out[["signif"]], Fun.Full.out[["signif"]]))

Procrustes <- read.csv("Outputs/Statistics/Procrustes_updated.csv") %>%
  select(-'X')

Procrustes_full2 <- rbind(Procrustes, Procrustes_full)

write.csv(Procrustes_full2, "Outputs/Statistics/Procrustes_fullrange.csv")

# Plot fit of procrustes analyses, Figure S1  
windows(width = 5, height = 4.5)
ProcrustesFit <- ggplot(Procrustes_full2, aes(x = factor(Spectral, levels = c('VNIR', 'SWIR', 'Full')), y = r, fill = factor(Field, levels =c("Taxonomic", "Functional")))) +
  geom_bar(stat = 'identity', position = position_dodge(), colour = 'black') +
  geom_text(aes(label = round(r, 2)), position = position_dodge(0.9), vjust = 1.5) +
  scale_fill_manual(values = c("Grey40", "Grey80")) +
  #geom_text(aes(label = r), position=position_dodge(width=0.9), vjust=-0.25, size=4)+
  scale_y_continuous(expand=c(0,0),breaks=seq(0, 1, by=0.25), limits=c(0,1))+ 
  scale_x_discrete(labels=c('VNIR', 'SWIR', 'Full Range'))+
  theme_classic() +
  ylab("Procrustean r")+
  xlab("")+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(face = 'bold', size = 10, colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 8, colour = 'black'),
        axis.title.y = element_text(face ='bold', size = 10),
        legend.text = element_text(face = 'bold', size = 10),
        legend.title = element_blank(),
        legend.position = c(0.85, 0.95)) 
ProcrustesFit

ggsave("Outputs/Figures/FigS1_ProcrustesFit_full.jpeg", width=5, height=4.5, dpi=600)

### 4.1.2__Pearson's Correlations ----------------------------------------------------
#   Examine if the individual PCs are associated with each other.
#   Visualize species, trait, wavelenght scores along PC axes

PCA_outputs_full2 <- PCA_outputs_full %>%
  pivot_wider(names_from = dimension, values_from = c(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8))

# Pearson's correlations of the first 3 PC axes
PCACor.Total <- cor(PCA_outputs_full2[,2:16]) # Correlate the PCs 1-3 of each dimension
colnames(PCACor.Total) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1", 'Full_PC1',
                            "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2", 'Full_PC2',
                            "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3", 'Full_PC3')
rownames(PCACor.Total) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1", 'Full_PC1',
                            "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2", 'Full_PC2',
                            "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3", 'Full_PC3')

PCACor.Test = cor.mtest(PCA_outputs_full2[,2:16], conf.level = 0.95) # Test correlations of the PCs 1-3 of each dimension
colnames(PCACor.Test[['p']]) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1", 'Full_PC1',
                                  "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2", 'Full_PC2',
                                  "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3", 'Full_PC3')
rownames(PCACor.Test[['p']]) <- c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1", 'Full_PC1',
                                  "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2", 'Full_PC2',
                                  "Tax_PC3", "Fun_PC3", "VNIR_PC3", "SWIR_PC3", 'Full_PC3')

# # Create a dataframe with correlation coefficients for the first 3 PCs across dimensions
PCACor.Total2 <- pull_lower_triangle(PCACor.Total) %>% 
  gather(Var2, r, 2:16) %>%
  rename(Var1 = rowname) %>%
  mutate(r = as.numeric(r)) %>%
  filter(!is.na(r)) %>%
  mutate(r = round(r, 2)) 

PCACor.Test2 <- pull_lower_triangle(PCACor.Test[["p"]]) %>% 
  gather(Var2, p, 2:16) %>%
  rename(Var1 = rowname) %>%
  mutate(p = as.numeric(p)) %>%
  filter(!is.na(p)) %>%
  mutate(p = round(p, 3)) 

PCACor.Total3 <- left_join(PCACor.Total2, PCACor.Test2) %>%
  filter(Var1 == "VNIR_PC1" | Var1 == "VNIR_PC2" | Var1 == "VNIR_PC3" |
           Var1 == "SWIR_PC1" | Var1 == "SWIR_PC2" | Var1 == "SWIR_PC3" |
           Var1 == "Full_PC1" | Var1 == "Full_PC2" | Var1 == "Full_PC3") %>%
  filter(Var2 == "Tax_PC1" | Var2 == "Tax_PC2" | Var2 == "Tax_PC3" |
           Var2 == "Fun_PC1" | Var2 == "Fun_PC2" | Var2 == "Fun_PC3" ) %>%
  rename(Spectra = Var1, Field = Var2)

write.csv(PCACor.Total3, "Outputs/Statistics/PCA_correlations_full.csv")

# Plot the Correlation matrix 
PCACor.Total4 <- left_join(PCACor.Total2, PCACor.Test2) %>%
  mutate(r2 = ifelse(p <= 0.05, r, NA))  # Remove r values for non-signif correlations -- FOR PLOTTING PURPOSES!!

PCACorTotalMatrix <- ggplot(data = PCACor.Total4, aes(Var2, Var1, fill = r))+
  scale_x_discrete(limits = c("Tax_PC1", "Fun_PC1", "VNIR_PC1", "SWIR_PC1", 'Full_PC1',
                              "Tax_PC2", "Fun_PC2", "VNIR_PC2", "SWIR_PC2", 'Full_PC2',
                              "Tax_PC3", "Fun_PC3", "VNIR_PC3", 'SWIR_PC3')) +
  scale_y_discrete(limits = c('Full_PC3',"SWIR_PC3","VNIR_PC3","Fun_PC3", "Tax_PC3",
                              'Full_PC2', "SWIR_PC2","VNIR_PC2","Fun_PC2", "Tax_PC2",
                              'Full_PC1', "SWIR_PC1","VNIR_PC1","Fun_PC1")) +
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

ggsave(plot = PCACorTotalMatrix, filename = 'Outputs/Figures/PCACor_Total_full.jpeg', width=6, height=6, dpi=300)


### 4.1.3__Understand PC axes ------------------------------------------------

Full.band.scores2 <- Full.band.scores %>%
  mutate(Region = ifelse(Wavelength <= 500, 'blue',
                         ifelse(Wavelength > 500 & Wavelength <= 600, 'green',
                                ifelse(Wavelength > 600 & Wavelength <= 730, 'red',
                                       ifelse(Wavelength > 730 & Wavelength <=780, 'rededge', 
                                              ifelse(Wavelength >780 & Wavelength <= 1400, 'NIR', 'SWIR'))))))

PC1.Full <- ggplot(Full.band.scores2, aes(x=Wavelength, y = PC1)) +
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
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'green'), mapping = aes(x=Wavelength, y = PC1), fill ='green') +
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'red'), mapping = aes(x=Wavelength, y = PC1), fill ='red') +
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'rededge'), mapping = aes(x=Wavelength, y = PC1), fill ='darkred') +
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'NIR'), mapping = aes(x=Wavelength, y = PC1), fill ='purple') +
  annotate('rect', xmin = 899, xmax = 957, ymin = -0.21, ymax = 0.21, colour = 'grey90', fill = 'grey90') +
  annotate('rect', xmin=c(1333,1790), xmax=c(1466,1950), ymin=c(-0.21,-0.21), ymax=c(0.21,0.21), colour = 'grey90', fill = 'grey90') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab("Wavelength (nm)") + 
  #  ylab("VNIR-Spectral PC1 \n (85%)") + 
  ylab('PC1 Scores') +
  labs(title = 'Full Range Spectral', subtitle = 'Variance Explained = 70%') +
  scale_y_continuous(limits = c(-0.21, 0.21), expand = c(0,0)) +
  scale_x_continuous(limits = c(400,2400), breaks = seq(400,2400,500)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(face='bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(size = 8, color='black'),
        axis.text.y = element_text(size = 8, color='black'),
        axis.title.y = element_text(face='bold', size = 10),
        axis.title.x = element_text(face='bold', size = 10),
        legend.title = element_text(face='bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5)
PC1.Full 

PC2.Full <- ggplot(Full.band.scores2, aes(x=Wavelength, y = PC2)) +
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
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'green'), mapping = aes(x=Wavelength, y = PC1), fill ='green') +
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'red'), mapping = aes(x=Wavelength, y = PC1), fill ='red') +
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'rededge'), mapping = aes(x=Wavelength, y = PC1), fill ='darkred') +
  # geom_area(data = Full.band.scores2 %>% filter(Region == 'NIR'), mapping = aes(x=Wavelength, y = PC1), fill ='purple') +
  annotate('rect', xmin = 899, xmax = 957, ymin = -0.21, ymax = 0.21, colour = 'grey90', fill = 'grey90') +
  annotate('rect', xmin=c(1333,1790), xmax=c(1466,1950), ymin=c(-0.21,-0.21), ymax=c(0.21,0.21), colour = 'grey90', fill = 'grey90') +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab("Wavelength (nm)") + 
  #  ylab("VNIR-Spectral PC1 \n (85%)") + 
  ylab('PC2 Scores') +
  labs(title = 'Full Range Spectral', subtitle = 'Variance Explained = 13%') +
  scale_y_continuous(limits = c(-0.21, 0.21), expand = c(0,0)) +
  scale_x_continuous(limits = c(400,2400), breaks = seq(400,2400,500)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(face='bold'),
        plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(size = 8, color = 'black'),
        axis.text.y = element_text(size = 8, color = 'black'),
        axis.title.y = element_text(face ='bold', size = 10),
        axis.title.x = element_text(face ='bold', size = 10),
        legend.title = element_text(face ='bold', size = 10),
        legend.position = 'bottom',
        legend.title.align=0.5)
PC2.Full 

windows(width = 7.08661, height = 6)
PCPlots <- ggarrange(PC1.VNIR, PC1.Full, PC1.SWIR, PC2.Full,
                     nrow = 2, ncol = 2,
                     common.legend = T, legend = 'bottom',
                     labels =  'AUTO')
PCPlots
ggsave(plot = PCPlots, filename = 'Outputs/Figures/FigS2_PC_plots3_full.jpeg', width=7.08661, height=6, dpi=600)


## 4.2__Full Spectral Diversity ---------------------------------------------

# This function standardizes the number of spectral points per VNIR and SWIR datasets, it is specific to the dataframes used here 
subsampleJoin <- function(VNIR, # Dataframe with VNIR spectral reflectance
                          SWIR, # Dataframe with SWIR spectral reflectance
                          min_points # Minimum number of spectral observations within area of interest (ie. plot)
                          ){
  VNIR2 <- VNIR %>%
    select(c(plot_field_id, '454.2': '970.77')) %>%
    group_by(plot_field_id) %>%
    sample_n({{min_points}}) %>%
    ungroup() %>%
    mutate(point_id = row_number(.), .before = plot_field_id)
  
  SWIR2 <- SWIR %>%
    select(c(plot_field_id, '972.5':'2412.5')) %>%
    group_by(plot_field_id) %>%
    sample_n({{min_points}}) %>%
    ungroup() %>%
    mutate(point_id = row_number(.), .before = plot_field_id)
  
  Full <- left_join(VNIR2, SWIR2) %>%
    select(-point_id)
  
  return(Full)
}

# Calculate coefficient of variation
CV_full <- replicate(999, # Number of resampling events, resampling occurs when fusing the spectra
                  subsampleJoin(Spec_VNIR, Spec_SWIR, 95) %>%  # Create fused spectra
                    calc.CV(plot_field_id, # Calculate diversity per plot
                            c('454.2':'2412.5'),  # Use full range of wavelengths to calulate diversity
                            FALSE, # Don't do rarefraction step, this is done when using spectra
                            NA)) # No rarefraction, so no internal resampling in diversity calculation.

# Create dataframe from list output
CV_full <- data.frame(CV_full[1,1], CV_full[2,]) 

# Average across each resampling event  
CV_full <- CV_full %>%  
  rowwise() %>%
  mutate(CV = mean(c_across(2:1000)), .after = plot_field_id) %>%
  select(c('plot_field_id', 'CV'))

# Calculate spectral variation
SV_full <- replicate(999,
                subsampleJoin(Spec_VNIR, Spec_SWIR, 95) %>% 
                  calc.SV(plot_field_id, # Calculate diversity per plot
                          c('454.2':'2412.5'),  # Use full range of wavelengths to calulate diversity
                          FALSE, # Don't do rarefraction step, this is done when using spectra
                          NA)) # No rarefraction, so no internal resampling in diversity calculation.

# Create dataframe from list output
SV_full <- data.frame(SV_full[1,1], SV_full[2,]) 

# Average across each resampling event
SV_full <- SV_full %>%  
  rowwise() %>%
  mutate(SV = mean(c_across(2:1000)), .after = plot_field_id) %>%
  select(c('plot_field_id', 'SV'))

# Add full range spectral diversity to diversity indices dataframe
# DivIndices <- read.csv(file = "Outputs/Data/DiversityIndices_updated.csv") %>%
#   select(-X)

DivIndices_full <- left_join(DivIndices, CV_full) %>%
  left_join(SV_full) %>%
  rename(SV.Full = SV,
         CV.Full = CV) %>%
  #select(-c('CHV.Casi', 'CHV.Sasi')) #%>%
  select(-c('CHV_VNIR', 'CHV_SWIR'))

write.csv(DivIndices_full, "Outputs/Data/DiversityIndices_full.csv", row.names = F)


### 4.2.1__Pearsons Correlation ----------------------------------------------
# Do correlations between all diversity metrics
DivCor.total <- cor(DivIndices_full[,-1], use="complete.obs") #Do correlation, where we exclude Na values (ie. plots where div index could be calculated)

# Name cols and rows of correlation matrix with appropriate diversity metric names
colnames(DivCor.total) <- c("q0", "q1", "q2", "J",
                            "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                            'CV_VNIR', 'SV_VNIR',"CV_SWIR",
                            'SV_SWIR', 'CV_Full',"SV_Full")
rownames(DivCor.total) <- c("q0", "q1", "q2", "J",
                            "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                            'CV_VNIR', 'SV_VNIR',"CV_SWIR",
                            'SV_SWIR', 'CV_Full',"SV_Full")

# Test the correlation relationships
DivCor.Test = cor.mtest(DivIndices_full[,-1], use="complete.obs", conf.level = 0.95)
colnames(DivCor.Test[['p']]) <- c("q0", "q1", "q2", "J",
                                  "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                                  'CV_VNIR', 'SV_VNIR',"CV_SWIR",
                                  'SV_SWIR', 'CV_Full',"SV_Full")
rownames(DivCor.Test[['p']]) <- c("q0", "q1", "q2", "J",
                                  "FRic", "FEve", "FDiv", "FDis", "RaoQ",
                                  'CV_VNIR', 'SV_VNIR',"CV_SWIR",
                                  'SV_SWIR', 'CV_Full',"SV_Full")

# Convert correlation matrix into dataframe
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
  filter(Var1 == "CV_VNIR" | Var1 == "CV_SWIR" |Var1 == "CV_Full" |
         Var1 == "SV_VNIR" | Var1 == "SV_SWIR" |Var1 == "SV_Full" ) %>%
  filter(Var2 == "q0" | Var2 == "q1" | Var2 == "q2" | Var2 == "J" |
         Var2 == "FRic" | Var2 == "FEve" | Var2 == "FDiv" | Var2 == "FDis" | Var2 == "RaoQ") %>%
  rename(Spectral = Var1, Field = Var2)

write.csv(DivCor, file = "Outputs/Statistics/DiversityCorrelations_full.csv")


# Summarize correlations by Sensor, Field Dimension, and Spectral Div metric and Plot
DivCor2 <- DivCor %>%
  separate(Spectral, c("Spec_metric","Spec_region"), sep = '[_]') %>%
  mutate(Field_dimension = ifelse(Field == "q0" |
                                  Field == "q1" | 
                                  Field == "q2" | 
                                  Field == "J", "Tax", "Func"))  %>%
  rename(Field_metric = Field)

## Plot relationships between spectral and field base diversity - Full is never better than just VNIR
# ggplot(DivCor, aes(x = Field, y = r, pch = Spectral, colour = Spectral)) +
#   geom_point(cex = 3) +
#   xlim('q0','q1','q2','J','FRic','FEve','FDiv','FDis','RaoQ')
# ggsave(filename = 'Outputs/Figures/DiversityCorrelations_metric_full.jpeg', width=5, height=4.5, dpi=600)

# Calculate average and sd, between field dimension and different spectral dimensions
DivCor2Sum <- DivCor2 %>%
  group_by(Spec_region, Field_dimension, Spec_metric) %>%
  mutate(mean = mean(r), se = sd(r)/sqrt(n()))

windows(width = 5, height = 4.5)
DivCorSum <- ggplot(DivCor2Sum, aes(x = factor(Field_dimension, level = c('Tax','Func')), y = mean, group = factor(Spec_metric, level = c('CV',"SV")))) +
  geom_bar(stat = 'identity', position = position_dodge(), aes(fill = factor(Spec_metric, level = c('CV', 'CHV', 'SV'))), colour = 'black') +
  #geom_text(aes(label = round(mean, 2)), position = position_dodge(0.9), vjust = 1.5) +
  geom_errorbar(aes(ymin = mean, ymax = mean +se), position = position_dodge(0.9), width = 0.2) +
  facet_wrap(~factor(Spec_region, levels = c('VNIR', 'SWIR', 'Full'), labels = c('VNIR','SWIR',"Full Range")), labeller = labeller(Sensor = c("Casi" = "VNIR", "Sasi" = "SWIR", "Full" = 'Full Range'))) +
  scale_fill_grey(start = 0.8, end = 0.2, name = 'Spectral Metric', labels = c('CV','SV')) +
  scale_x_discrete("", labels = c("Tax", 'Fun')) +
  scale_y_continuous("Mean Correlation with Spectral Diversity \n(Pearson's r)", expand = c(0,0), limits = c(-0,1), breaks = c(0,.25, .5,.75,1))+
  theme_classic() +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(face ='bold', size=10, color='black'),
        axis.text.y = element_text(size = 8, colour = 'black'),
        axis.title.y = element_text(face = 'bold', size = 10),
        legend.title = element_text(face = 'bold', size = 10),
        legend.text.align = 0.5,
        legend.position = 'bottom',
        #legend.box.background = element_rect(colour = "black"),
        legend.title.align = 0.5,
        strip.text = element_text(face = 'bold', size = 12))
DivCorSum

ggsave(plot = DivCorSum, filename = 'Outputs/Figures/FigS3_DiversityCorrelations_Full.jpeg', width=5, height=4.5, dpi=600)


# Create wide format
DivCor3 <- DivCor2 %>%
  mutate(r = round(r, 3)) %>%
  mutate(p = round(p, 3)) %>%
  pivot_wider(names_from = Spec_metric, values_from = c("r", "p"))

write.csv(DivCor3, "Outputs/Statistics/DiversityCorrelations_wide_Full.csv")

# Calculate differences between correlation with full vs spectral regions, Table S
DivCor4 <- DivCor3.2 %>%  
  dplyr::select(-c('p_CV', 'p_SV')) %>%
  pivot_wider(names_from = Spec_region, values_from = c('r_CV', 'r_SV')) %>%
  group_by(Field_metric) %>%
  mutate(dif_CV_VNIR = sum(r_CV_Full - r_CV_VNIR),
         dif_CV_SWIR = sum(r_CV_Full - r_CV_SWIR),
         dif_SV_VNIR = sum(r_SV_Full - r_SV_VNIR),
         dif_SV_SWIR = sum(r_SV_Full - r_SV_SWIR)) %>%
  ungroup() %>%
  dplyr::select(-c('r_CV_VNIR','r_CV_SWIR','r_SV_VNIR','r_SV_SWIR')) 

write.csv(DivCor4, "Outputs/Statistics/TableS2_DiversityCorrelations_wide_Full.csv")

