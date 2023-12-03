# Calculation of fecal sample alpha diversity for "Predictors and 
# consequences of diet composition in a declining generalist aerial insectivore"

# Written by: Jenny Uehling and Conor Taff
# Last updated: 11/5/2023
# Run under R Studio with R version 4.3.1 on a Mac OS

# This code takes the phyloseq object coi_ps2, created in data_filtering_script.R,
# and performs calculations to determine the alpha diversity of each fecal sample.

# Load libraries ---------------------------------------------------------------

pacman::p_load("tidyverse", "plyr", "dplyr", "phyloseq", "vegan", "here",
               "data.table")

# tidyverse, plyr, and dplyr for data wrangling
# phyloseq & vegan for community analyses and plotting
# here for file reference by relative paths
# data.table for creating data tables

# Will need to install these libraries first if not yet installed
# To install phyloseq, do the following:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

# Load phyloseq object ---------------------------------------------------------
coi_ps2 <- readRDS("2_modified_data/coi_ps2.rds")

# Load s_info ------------------------------------------------------------------
# This file contains information about each fecal sample that is needed for future
# analyses.
s_info <- read.csv("2_modified_data/s_info.csv")

# Agglomerate taxa -------------------------------------------------------------

# Agglomerate to family
coi_fam2 <- tax_glom(coi_ps2, taxrank = "family")

# Filter samples, ASVs based on read counts ------------------------------------

# Create a record of the sequencing depth of each sample before pruning
depth_preprune <- data.frame(as(sample_data(coi_fam2), "data.frame"),
                                     TotalReads = sample_sums(coi_fam2), keep.rownames = TRUE)
depth_preprune <- subset(depth_preprune, select = c(sampleID, Individual_Band, Age, Site, Nest, TotalReads)) # for readability

# Remove ASVs with less than 10 reads across all samples
coi_ps2 <- prune_taxa(taxa_sums(coi_fam2) > 10, coi_fam2) # Based on Hoenig et al. 2020 and Forsman et al. 2022

# Remove samples with less than 100 total reads per sample
coi_ps2 <- prune_samples(sample_sums(coi_ps2) >= 100, coi_ps2)

# Create a record of the sequencing depth of each sample after filtering
depth_postprune <- data.frame(as(sample_data(coi_ps2), "data.frame"),
                              TotalReads = sample_sums(coi_ps2), keep.rownames = TRUE)
depth_postprune <- subset(depth_postprune, select = c(sampleID, Individual_Band, Age, Site, Nest, TotalReads)) # for readability

# Explore alpha diversity ------------------------------------------------------

# We will use coi_ps2 for alpha diversity calculations

## Calculate diversity indices without rarefying data ------------------------

# Put data in correct form for vegan package
# Rows are sampleIDs, columns are families, each cell is number of reads
coi_ps2_psmelt <- psmelt(coi_ps2)
coi_ps2_psmelt <- subset(coi_ps2_psmelt, select = c(sampleID, family, Abundance))
coi_ps2_psmelt_wide <- coi_ps2_psmelt %>% pivot_wider(
  names_from = family,
  values_from = Abundance
)
coi_ps2_psmelt_wide <- coi_ps2_psmelt_wide %>% remove_rownames %>% column_to_rownames(var="sampleID") # Make row names sampleIDs

# Calculate Shannon index
shannon_psmelt <- diversity(coi_ps2_psmelt_wide, index = "shannon")
shannon_psmelt <- data.frame(shannon_psmelt)
shannon_psmelt <- rownames_to_column(shannon_psmelt, "sampleID")

# Calculate Simpson index
simpson_psmelt <- diversity(coi_ps2_psmelt_wide, index = "simpson")
simpson_psmelt <- data.frame(simpson_psmelt)
simpson_psmelt <- rownames_to_column(simpson_psmelt, "sampleID")

# Calculate number of families
num_fam_psmelt <- data.frame(matrix(NA, ncol = 2, nrow = length(unique(coi_ps2_psmelt$sampleID))))
sampleID_list <- unique(coi_ps2_psmelt$sampleID)
for (i in 1:length(sampleID_list)){
  sam <- sampleID_list[i]
  sub <- coi_ps2_psmelt[coi_ps2_psmelt$sampleID == sam ,]
  num_fam <- nrow(sub[sub$Abundance > 0 ,])
  num_fam_psmelt[i,1] <- sam
  num_fam_psmelt[i,2] <- num_fam
}
num_fam_psmelt <- dplyr::rename(num_fam_psmelt,
                                sampleID = X1, num_fam_psmelt = X2)

## Rarefy samples ------------------------------------------------------------

# Rarefy to even depth of 150, 100, and 50
coi_ps2_rar150_raw <- rarefy_even_depth(coi_ps2, sample.size = 150, rngseed = 92)
coi_ps2_rar100_raw <- rarefy_even_depth(coi_ps2, sample.size = 100, rngseed = 92)
coi_ps2_rar50_raw <- rarefy_even_depth(coi_ps2, sample.size = 50, rngseed = 92)

## Calculate diversity indices from rar150 -----------------------------------

# Put data in correct form for vegan package
# Rows are sampleIDs, columns are families, each cell is number of reads
coi_ps2_rar150 <- psmelt(coi_ps2_rar150_raw)
coi_ps2_rar150 <- subset(coi_ps2_rar150, select = c(sampleID, family, Abundance))
coi_ps2_rar150_wide <- coi_ps2_rar150 %>% pivot_wider(
  names_from = family,
  values_from = Abundance
)
coi_ps2_rar150_wide <- coi_ps2_rar150_wide %>% remove_rownames %>% column_to_rownames(var="sampleID") # Make row names sampleIDs

# Calculate Shannon index
shannon_150 <- diversity(coi_ps2_rar150_wide, index = "shannon")
shannon_150 <- data.frame(shannon_150)
shannon_150 <- rownames_to_column(shannon_150, "sampleID")

# Calculate Simpson index
simpson_150 <- diversity(coi_ps2_rar150_wide, index = "simpson")
simpson_150 <- data.frame(simpson_150)
simpson_150 <- rownames_to_column(simpson_150, "sampleID")

# Calculate number of families
num_fam_150 <- data.frame(matrix(NA, ncol = 2, nrow = length(unique(coi_ps2_rar150$sampleID))))
sampleID_list <- unique(coi_ps2_rar150$sampleID)
for (i in 1:length(sampleID_list)){
  sam <- sampleID_list[i]
  sub <- coi_ps2_rar150[coi_ps2_rar150$sampleID == sam ,]
  num_fam <- nrow(sub[sub$Abundance > 0 ,])
  num_fam_150[i,1] <- sam
  num_fam_150[i,2] <- num_fam
}
num_fam_150 <- dplyr::rename(num_fam_150,
                             sampleID = X1, num_fam_150 = X2)

## Calculate diversity indices from rar100 -----------------------------------

# Put data in correct form for vegan package
# Rows are sampleIDs, columns are families, each cell is number of reads
coi_ps2_rar100 <- psmelt(coi_ps2_rar100_raw)
coi_ps2_rar100 <- subset(coi_ps2_rar100, select = c(sampleID, family, Abundance))
coi_ps2_rar100_wide <- coi_ps2_rar100 %>% pivot_wider(
  names_from = family,
  values_from = Abundance
)
coi_ps2_rar100_wide <- coi_ps2_rar100_wide %>% remove_rownames %>% column_to_rownames(var="sampleID") # Make row names sampleIDs

# Calculate Shannon index
shannon_100 <- diversity(coi_ps2_rar100_wide, index = "shannon")
shannon_100 <- data.frame(shannon_100)
shannon_100 <- rownames_to_column(shannon_100, "sampleID")

# Calculate Simpson index
simpson_100 <- diversity(coi_ps2_rar100_wide, index = "simpson")
simpson_100 <- data.frame(simpson_100)
simpson_100 <- rownames_to_column(simpson_100, "sampleID")

# Calculate number of families
num_fam_100 <- data.frame(matrix(NA, ncol = 2, nrow = length(unique(coi_ps2_rar100$sampleID))))
sampleID_list <- unique(coi_ps2_rar100$sampleID)
for (i in 1:length(sampleID_list)){
  sam <- sampleID_list[i]
  sub <- coi_ps2_rar100[coi_ps2_rar100$sampleID == sam ,]
  num_fam <- nrow(sub[sub$Abundance > 0 ,])
  num_fam_100[i,1] <- sam
  num_fam_100[i,2] <- num_fam
}
num_fam_100 <- dplyr::rename(num_fam_100,
                             sampleID = X1, num_fam_100 = X2)

## Calculate diversity indices from rar50 ------------------------------------

# Put data in correct form for vegan package
# Rows are sampleIDs, columns are families, each cell is number of reads
coi_ps2_rar50 <- psmelt(coi_ps2_rar50_raw)
coi_ps2_rar50 <- subset(coi_ps2_rar50, select = c(sampleID, family, Abundance))
coi_ps2_rar50_wide <- coi_ps2_rar50 %>% pivot_wider(
  names_from = family,
  values_from = Abundance
)
coi_ps2_rar50_wide <- coi_ps2_rar50_wide %>% remove_rownames %>% column_to_rownames(var="sampleID") # Make row names sampleIDs

# Calculate Shannon index
shannon_50 <- diversity(coi_ps2_rar50_wide, index = "shannon")
shannon_50 <- data.frame(shannon_50)
shannon_50 <- rownames_to_column(shannon_50, "sampleID")

# Calculate Simpson index
simpson_50 <- diversity(coi_ps2_rar50_wide, index = "simpson")
simpson_50 <- data.frame(simpson_50)
simpson_50 <- rownames_to_column(simpson_50, "sampleID")

# Calculate number of families
num_fam_50 <- data.frame(matrix(NA, ncol = 2, nrow = length(unique(coi_ps2_rar50$sampleID))))
sampleID_list <- unique(coi_ps2_rar50$sampleID)
for (i in 1:length(sampleID_list)){
  sam <- sampleID_list[i]
  sub <- coi_ps2_rar50[coi_ps2_rar50$sampleID == sam ,]
  num_fam <- nrow(sub[sub$Abundance > 0 ,])
  num_fam_50[i,1] <- sam
  num_fam_50[i,2] <- num_fam
}
num_fam_50 <- dplyr::rename(num_fam_50,
                            sampleID = X1, num_fam_50 = X2)

## Correlation between diversity metrics -------------------------------------

# Create dataframe with all metrics calculated across different levels of rarefaction
data_list <- list(shannon_psmelt, simpson_psmelt, num_fam_psmelt, 
                  shannon_150, simpson_150, num_fam_150,
                  shannon_100, simpson_100, num_fam_100,
                  shannon_50, simpson_50, num_fam_50)

alpha_div <- data_list %>% reduce(full_join, by = "sampleID")

# Shannon
mod <- lm(shannon_psmelt~shannon_150, data = alpha_div)
summary(mod)
mod <- lm(shannon_psmelt~shannon_100, data = alpha_div)
summary(mod)
mod <- lm(shannon_psmelt~shannon_50, data = alpha_div)
summary(mod)

# Simpson
mod <- lm(simpson_psmelt~simpson_150, data = alpha_div)
summary(mod)
mod <- lm(simpson_psmelt~simpson_100, data = alpha_div)
summary(mod)
mod <- lm(simpson_psmelt~simpson_50, data = alpha_div)
summary(mod)

# Number of families
mod <- lm(num_fam_psmelt~num_fam_150, data = alpha_div)
summary(mod)
mod <- lm(num_fam_psmelt~num_fam_100, data = alpha_div)
summary(mod)
mod <- lm(num_fam_psmelt~num_fam_50, data = alpha_div)
summary(mod)

# Simpson and Shannon
mod <- lm(simpson_psmelt~shannon_psmelt, data = alpha_div)
summary(mod)

# It doesn't seem like rarefying does much to change Shannon and Simpson,
# as the ajdusted R-squared values are very high. There are big differences for
# richness (number of families), but this is probably not as useful of a metric.

# Simpson and Shannon are also very closely related (adjusted R-squared: 0.9344).
# Use Simpson for future analyses. Save dataframe just with Simpson.

alpha_div <- subset(alpha_div, select = c(sampleID, simpson_psmelt))

# Write csv to use in data analysis script ------------------------------------

write.csv(alpha_div, here("2_modified_data/alpha_div.csv"))
