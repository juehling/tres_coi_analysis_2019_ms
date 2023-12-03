# Calculation of proportion of swallow diets composed of aquatic insects for "Predictors and 
# consequences of diet composition in a declining generalist aerial insectivore"

# Written by: Jenny Uehling and Conor Taff
# Last updated: 12/3/2023
# Run under R Studio with R version 4.3.1 on a Mac OS

# This code takes the phyloseq object coi_ps2, created in data_filtering_script.R,
# and performs calculations to determine the percentage of each fecal sample
# composed of aquatic insects.

# Load libraries ---------------------------------------------------------------

pacman::p_load("tidyverse", "plyr", "dplyr", "phyloseq", "vegan", "here",
               "data.table", "ggpubr", "grid")

# tidyverse, plyr, and dplyr for data wrangling
# phyloseq & vegan for community analyses and plotting
# here for file reference by relative paths
# data.table for creating data tables
# ggpubr and grid for plotting

# Will need to install these libraries first if not yet installed
# To install phyloseq, do the following:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("phyloseq")

# Load phyloseq object ---------------------------------------------------------

coi_ps2 <- readRDS("2_modified_data/coi_ps2.rds")

# Load s_info ------------------------------------------------------------------

# This file contains information about each fecal sample that is needed for future
# analyses. It was created in data_filtering_script.R
s_info <- read.csv("2_modified_data/s_info.csv")

# Add in family life history information ---------------------------------------

# Import the rest of the information we have about the larval stage of each
# insect at the family level.

# Merge family to life history
life_history_family <- read.csv(here("4_other_outputs", "unique_families_aquatic_terrestrial_IDs.csv")) # prepared outside of R

# Join tax_table from coi_ps2 and family life history .csv
otu_lh <- plyr::join(as.data.frame(tax_table(coi_ps2)), life_history_family, by = "family", type = "left")

# For readability, subset just to columns needed in analysis
otu_lh <- subset(otu_lh, select = c(kingdom, phylum, class, order, family, genus, species, life_history))

rownames(otu_lh) <- rownames(tax_table(coi_ps2))

# When agglomerating to a specific taxonomic level, phyloseq deletes everything
# after that taxonomic level in the order of the columns. Therefore, we have to put
# "life_history" before "family."
col_order <- c("kingdom", "phylum", "class", "order", "life_history", "family", "genus", "species")
otu_lh <- otu_lh[, col_order]

# Save all life history information in new phyloseq object ---------------------

# Import otu_lh as new tax_table in coi_ps2
coi_ps2 <- phyloseq(
  otu_table(coi_ps2),
  tax_table(as.matrix(otu_lh)),
  sample_data(coi_ps2)
)

# Agglomerate taxa -------------------------------------------------------------

# Agglomerate to family
coi_fam2 <- tax_glom(coi_ps2, taxrank = "family")

# Rarefaction plots ------------------------------------------------------------

# Plot rarefaction curves with rarecurve function from vegan package.
# We're doing this to figure out a reasonable number of minimum reads for each sample.
# We'll remove any samples with less than a certain number of reads.
# Split out into nestlings and adults for readability

# Create curves with coi_fam2 object

# Nestling curves
coi_fam2_nestling <- subset_samples(coi_fam2, Adult_or_Nestling == "Nestling")
# Must first put OTUs in columns as per vegan standards.
veg_tab_nestling <- otu_table(coi_fam2_nestling)
class(veg_tab_nestling) <- "matrix"
rarecurve(t(veg_tab_nestling),step=10,label=FALSE,cex=0.5,xlim=c(0,500), ylab="Families")
title("Nestlings")

# Adult curves
coi_fam2_adult <- subset_samples(coi_fam2, Adult_or_Nestling == "Adult")
# Must first put OTUs in columns as per vegan standards.
veg_tab_adult <- otu_table(coi_fam2_adult)
class(veg_tab_adult) <- "matrix"
rarecurve(t(veg_tab_adult),step=10,label=FALSE,cex=0.5,xlim=c(0,500), ylab="Families")
title("Adults")

# Filter samples, ASVs based on read counts, relative abundance ----------------

# Create a record of the sequencing depth of each sample before pruning
depth_preprune <- data.frame(as(sample_data(coi_fam2), "data.frame"),
                                     TotalReads = sample_sums(coi_fam2), keep.rownames = TRUE)
depth_preprune <- subset(depth_preprune, select = c(sampleID, Individual_Band, Age, Site, Nest, TotalReads)) # for readability

# Remove ASVs with 10 or fewer reads across all samples
coi_ps2 <- prune_taxa(taxa_sums(coi_fam2) > 10, coi_fam2) # Based on Hoenig et al. 2020 and Forsman et al. 2022

# Remove samples with less than 100 total reads
coi_ps2 <- prune_samples(sample_sums(coi_ps2) >= 100, coi_ps2)

# Create a record of the sequencing depth of each sample before transforming to relative abundance
depth_postprune <- data.frame(as(sample_data(coi_ps2), "data.frame"),
                              TotalReads = sample_sums(coi_ps2), keep.rownames = TRUE)
depth_postprune <- subset(depth_postprune, select = c(sampleID, Individual_Band, Age, Site, Nest, TotalReads)) # for readability

# Transform to relative abundance
coi_ra2 <- transform_sample_counts(coi_ps2, function(x) x / sum(x))

# Transform to occurrence
coi_occ <- transform_sample_counts(coi_ra2, function(x) ceiling(x))

# Limit to families in 20% of samples
coi_20 <- prune_taxa(genefilter_sample(coi_occ, filterfun_sample(function(x) x > 0), A = 0.2 * nsamples(coi_occ)), coi_occ)

# Examine sequencing depth differences pre vs. post filtering ------------------

# Here, we examine how the ratio of arthropods to other living things compared in
# our sequences, and how filtering changed our sequencing depths. We use .csv 
# files created and saved in the data_filtering_script.R, as well as dataframes
# created earlier in this script.

  ## Arthropod reads vs. all reads comparison ----------------------------------

  # Read in sequencing depths without arthropod filter, created in data_filtering_script.R
  depth_no_arth_filter <- read.csv("2_modified_data/depth_no_arth_filter.csv")

  # Rename column to indicate all reads (including arthropods)
  depth_no_arth_filter <- dplyr::rename(depth_no_arth_filter,
                                      TotalReadsAll = TotalReads)
  
  # Call day 6 nestling Individual_Band values "NA" so that the dataframes merge
  # properly (they are NA in depth_no_arth_filter)
  depth_preprune$Individual_Band[depth_preprune$Age == "6"] <- NA
  
  # Merge all reads dataset with only arthropods dataset
  depth_examine <- merge(depth_preprune, depth_no_arth_filter)
  
  # Examine the percent of reads retained when filtering just to arthropod
  depth_examine$percent_retained <- depth_examine$TotalReads/depth_examine$TotalReadsAll
  mean(depth_examine$percent_retained)
  range(depth_examine$percent_retained)
  
  # Sum up total number of reads after filtering to arthropods
  sum(depth_preprune$TotalReads)
  
  # Determine the range of reads across samples
  range(depth_preprune$TotalReads)
  
  # Determine mean and median number of reads per sample
  mean(depth_preprune$TotalReads)
  median(depth_preprune$TotalReads)
  
  # Plot the total reads without filtering to arthropod vs. with filtering to arthropod
  p <- ggplot(depth_examine, aes(x=TotalReadsAll, y = TotalReads)) + geom_point() +
    theme_classic() + xlab("Total reads without filtering to Arthropod") +
    ylab("Total reads filtering to Arthropod") +
    theme(axis.title=element_text(size = 18))
  p

  ## Examine what proportion of each sample was Tachycineta reads ---------------------
  
  # Tachycineta is the genus for tree swallows. How many of our sequences were
  # from the animal from which the feces were produced?
  
  # Read in phyloseq object from before filtering to arthropod
  coi_ps <- readRDS("2_modified_data/coi_ps.rds")

  # Pull out just reads with genus Tachycineta
  coi_tres <- subset_taxa(coi_ps, genus == "Tachycineta")
  
  # Save sequencing depth
  depth_tres <- data.frame(as(sample_data(coi_tres), "data.frame"),
                         TotalReads = sample_sums(coi_tres), keep.rownames = TRUE)
  depth_tres <- subset(depth_tres, select = c(sampleID, TotalReads)) # for readability and ease of merging
  
  # Rename column so it can be merged with other datasets
  depth_tres <- dplyr::rename(depth_tres,
                            TotalReadsTRES = TotalReads)
  
  # Merge with other sequencing depth datasets
  depth_examine_tres <- merge(depth_examine, depth_tres, all.x = TRUE, all.y = FALSE)
  
  # Examine percent of each sample composed of Tachycineta reads
  depth_examine_tres$percent_tres <- depth_examine_tres$TotalReadsTRES/depth_examine_tres$TotalReadsAll
  
  # Determine mean percentage of Tachycineta reads
  mean(depth_examine_tres$percent_tres)

# Convert phyloseq objects to data frames  -------------------------------------

# At this point, we've done all of the filtering/organizing/plotting/modeling that
# we want with phyloseq.

# We will now convert the phyloseq objects into data frames so that we can more
# easily work with them.

# psmelt makes a phyloseq object into a data frame

plot_ra <- psmelt(coi_ra2)

plot_occ <- psmelt(coi_occ)

coi_20_df <- psmelt(coi_20)

# Identify families' mean relative abundance, other stats about dataset --------

# Calculate number of samples
num_sam <- length(unique(plot_ra$Sample))

# Extract raw families
family_raw <- unique(plot_ra$family)
family_raw <- as.data.frame(family_raw)
family_raw <- dplyr::rename(family_raw,
                         family = family_raw)

# Add in life history information
family_raw <- merge(family_raw, life_history_family, by = "family", all.x = TRUE, all.y = FALSE)
table(family_raw$life_history)

# Calculate percent of families with "both" classification
17/171

# Create data frame to store mean relative abundance information
mean_rel_ab_families <- data.frame(matrix(NA, ncol = 2, nrow = length(family_raw$family)))

for (i in 1:length(family_raw$family)){
  fam <- family_raw$family[i]
  plot_ra_fam <- plot_ra[plot_ra$family == fam ,]
  tot <- sum(plot_ra_fam$Abundance)
  mean_relab <- tot/num_sam
  mean_rel_ab_families[i,1] <- fam
  mean_rel_ab_families[i,2] <- mean_relab
}
# Examine mean_rel_ab_families table to determine mean values for each family.

# Summarize number of samples across categories in final dataset ---------------

# The information extracted in this section is reported in the text of the manuscript
# or was used to populate Table S1 and Table S2

sampleID <- unique(plot_ra$sampleID)
samples <- as.data.frame(sampleID)
samples_info <- merge(samples, s_info)
uniq_indiv <- as.data.frame(unique(plot_ra$Individual_Band))
table(samples_info$Adult_or_Nestling)
table(samples_info$Age)
table(samples_info$Sex)

# Break down the number of samples per individual by adult and nestling
# Adults
adults_samples_info <- samples_info[samples_info$Adult_or_Nestling == "Adult" ,]
table(adults_samples_info$Sex, adults_samples_info$Site, adults_samples_info$Capture_Number)
table(adults_samples_info$Capture_Number, adults_samples_info$Sex)
indiv_samples_adults <- as.data.frame(table(adults_samples_info$Individual_Band))
table(indiv_samples_adults$Freq)
# Nestlings
nestlings_samples_info <- samples_info[samples_info$Adult_or_Nestling == "Nestling" ,]
table(nestlings_samples_info$Age)
table(nestlings_samples_info$Age, nestlings_samples_info$Site)
# Figure out how many samples were taken per box
nest_boxes <- nestlings_samples_info[c("Site", "site_box_year")]
nest_boxes_unique <- unique(nest_boxes)
table(nest_boxes_unique$Site)
indiv_samples_nestlings <- as.data.frame(table(nestlings_samples_info$Individual_Band))
table(indiv_samples_nestlings$Freq) # 151 + 25 = 176 banded nestlings with samples

# Plot Patterns: Taxonomic Groups ----------------------------------------------

## Examine occurrence of families found in 20% or more of samples split by age ----- 

# Figure out how many families occurred in 20% or more of the samples
table(coi_20_df$sampleID)

# Calculate number of nestlings and adults to use for determining percentages of samples
coi_20_df_uniq <- dplyr::select(coi_20_df, c("sampleID", "Adult_or_Nestling"))
coi_20_df_uniq <- unique(coi_20_df_uniq)
num_ad <- nrow(coi_20_df_uniq[coi_20_df_uniq$Adult_or_Nestling == "Adult" ,])
num_nestl <- nrow(coi_20_df_uniq[coi_20_df_uniq$Adult_or_Nestling == "Nestling" ,])

# Plot adults and nestlings together
p_all <- ggplot(coi_20_df, aes(x=family, y=(Abundance/(num_ad + num_nestl)), fill = life_history)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  ylab("Percentage of samples with arthropod family (%)") +
  scale_fill_manual(values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4")) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  theme(plot.title = element_text(size = 14))

# Figure out order for families on the x axis
fams_20 <- unique(coi_20_df$family)
order_for_axis <- data.frame(matrix(NA, ncol = 2, nrow = length(fams_20)))
num_sam <- length(coi_20_df_uniq$sampleID)

for (i in 1:length(fams_20)){
  fam <- fams_20[i]
  coi_20_df_fam <- coi_20_df[coi_20_df$family == fam ,]
  tot <- sum(coi_20_df_fam$Abundance)
  percent_samples <- tot/num_sam
  order_for_axis[i,1] <- fam
  order_for_axis[i,2] <- percent_samples
}

# Order dataframe based on relative abundance of families
order_for_axis <- order_for_axis[order(-order_for_axis$X2) ,]

# Save for use in x-axis ordering
order_for_axis_list <- as.character(order_for_axis$X1)

# Plot families found in over 20% of samples, in adults and nestlings
p_adults <- ggplot(coi_20_df[coi_20_df$Adult_or_Nestling == "Adult" ,], aes(x=family, y=(Abundance/num_ad), fill = life_history)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_x_discrete(limits = order_for_axis_list) + 
  ylab("Proportion of samples with family") +
  scale_fill_manual(values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4")) +
  theme(legend.title = element_blank()) +
  ggtitle("Adults") +
  theme(plot.title = element_text(size = 14))

p_nestlings <- ggplot(coi_20_df[coi_20_df$Adult_or_Nestling == "Nestling" ,], aes(x=family, y=(Abundance/num_nestl), fill = life_history)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_x_discrete(limits = order_for_axis_list) + 
  xlab("Family") +
  ylab("Proportion of samples with family") +
  scale_fill_manual(values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4")) +
  ggtitle("Nestlings") +
  theme(plot.title = element_text(size = 14))

# Combine figures
fig_a <- ggarrange(p_adults + rremove("ylab") + rremove("xlab") +
                   theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         plot.margin = margin(b = 1, l = 7)),
                 p_nestlings + rremove("ylab") + rremove("xlab") +
                   theme(plot.margin = margin(t = 0.001, l = 7)),
                 ncol = 1, 
                 common.legend = TRUE, legend = "right",
                 align = "h")

fig_a <- annotate_figure(fig_a, left = textGrob("Proportion of samples with family", rot = 90, vjust = 1, gp = gpar(cex = 2)),
                bottom = textGrob("Arthropod family", gp = gpar(cex = 2)))

fig_a

# Save combined figure
ggsave(here("3_r_scripts/figs/figs_descriptive/common_families.pdf"), width = 7, height = 7, device = "pdf")

## Examine relative abundance of top 15 families with highest relative abundance split by age -----

# Order dataframe based on relative abundance of families
order_for_axis <- mean_rel_ab_families[order(-mean_rel_ab_families$X2) ,]

# Cut off only 15 most common families
order_for_axis <- order_for_axis[1:15,]

# Save for use in x-axis ordering
order_for_axis_list <- as.character(order_for_axis$X1)

# Calculate number of nestlings and adults to use for determining relative abundance
plot_ra_uniq <- dplyr::select(plot_ra, c("sampleID", "Adult_or_Nestling"))
plot_ra_uniq <- unique(plot_ra_uniq)
num_ad <- nrow(plot_ra_uniq[plot_ra_uniq$Adult_or_Nestling == "Adult" ,])
num_nestl <- nrow(plot_ra_uniq[plot_ra_uniq$Adult_or_Nestling == "Nestling" ,])

# Plot families with highest relative abundance, in adults and nestlings
p_adults <- ggplot(plot_ra[plot_ra$Adult_or_Nestling == "Adult" ,], aes(x=family, y=(Abundance/num_ad), fill = life_history)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_x_discrete(limits = order_for_axis_list) +
  scale_y_continuous(limits = c(0, 0.30)) +
  ylab("Average relative read abundance") +
  scale_fill_manual(values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4")) +
  theme(legend.title = element_blank()) +
  ggtitle("Adults") +
  theme(plot.title = element_text(size = 14))
# Note: this displays a warning message because we're only showing the most common families, and the rest are removed.

p_nestlings <- ggplot(plot_ra[plot_ra$Adult_or_Nestling == "Nestling" ,], aes(x=family, y=(Abundance/num_nestl), fill = life_history)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 12)) +
  scale_x_discrete(limits = order_for_axis_list) +
  scale_y_continuous(limits = c(0, 0.30)) +
  xlab("Family") +
  ylab("Average relative abundance") +
  scale_fill_manual(values = c("aquatic" = "skyblue1", "terrestrial" = "palegreen3", "both" = "tan4")) +
  theme(legend.title = element_blank()) +
  ggtitle("Nestlings") +
  theme(plot.title = element_text(size = 14))
# Note: this displays a warning message because we're only showing the most common families, and the rest are removed.

# Combine figures
fig_b <- ggarrange(p_adults + rremove("ylab") + rremove("xlab") +
                   theme(axis.text.x = element_blank(),
                         axis.ticks.x = element_blank(),
                         axis.title.x = element_blank(),
                         plot.margin = margin(b = 1, l = 7)),
                 p_nestlings + rremove("ylab") + rremove("xlab") +
                   theme(plot.margin = margin(t = 0.001, l = 7)),
                 ncol = 1, 
                 common.legend = TRUE, legend = "right",
                 align = "h")

fig_b <- annotate_figure(fig_b, left = textGrob("Average relative read abundance of family", rot = 90, vjust = 1, gp = gpar(cex = 2)),
                bottom = textGrob("Arthropod family", gp = gpar(cex = 2)))

fig_b

# Save combined figures
ggsave(here("3_r_scripts/figs/figs_descriptive/common_families_rel_ab.pdf"), width = 5, height = 7, device = "pdf")

## Combine figures for families found in over 20% of samples and families with highest relative abundance ----- 

# Opened "common_families.pdf" and "common_families_rel_ab.pdf" in Adobe Illustrator.

# Combined figures into single PDF and made the following modifications:

  # - Deleted duplicate legend
  # - Increased size of text in legend, made legend central
  # - Deleted duplicate x-axis label, centered label ("Arthropod family") and made it larger
  # - Moved y-axis labels to be more centered
  # - Ensured that axes in Figures A and B were at the same levels.
  # - Added "A" and "B" labels

# Exporting as TIF file results in best visuals.

# Calculate percent aquatic & number of families in each sample ----------------

## Percent aquatic using relative abundance ----------------------------------
# For this first calculation, we will also save the number of families per sample

# Identify the unique samples
sample <- unique(plot_ra$sampleID)

# Create an empty data frame to store relative abundance and number of families information
aquatic <- data.frame(matrix(NA, ncol = 3, nrow = length(sample)))

for (i in 1:length(sample)){
  sam <- sample[i] # identify sample
  list <- plot_ra[plot_ra$sampleID == sam ,] # pull out all records from plot_ra that have that sample ID
  unique_fam <- unique(list$family[list$Abundance != "0"])
  num_fam <- length(unique_fam)
  list_aq <- list[list$life_history == "aquatic" ,] # of those records, pull out only those that have aquatic life histories
  rel_ab <- sum(list_aq$Abundance) # sum up all of the percentages of aquatic insects for that sample
  aquatic[i,1] <- sam # save the sample name
  aquatic[i,2] <- rel_ab # save the total aquatic relative abundance
  aquatic[i,3] <- num_fam # save the number of families in the sample
}

# Name columns
names(aquatic)[names(aquatic) == "X1"] <- "sampleID"
names(aquatic)[names(aquatic) == "X2"] <- "percent_aquatic_ra"
names(aquatic)[names(aquatic) == "X3"] <- "num_fam"

# Add in sample information
aquatic <- merge(aquatic, s_info, by = "sampleID", all.x = TRUE, all.y = FALSE)

## Percent aquatic using occurrence --------------------------------------------

# Identify the unique samples
sample <- unique(plot_occ$sampleID)

# Create an empty data frame to store occurrence information
occ_aq <- data.frame(matrix(NA, ncol = 2, nrow = length(sample)))

for (i in 1:length(sample)){
  sam <- sample[i] # identify sample
  list <- plot_occ[plot_occ$sampleID == sam ,] # pull out all records from plot_occ that have that sample ID
  list <- list[list$Abundance != 0 ,] # Only include families in this calculation if they're actually present in the sample
  list_aq <- list[list$life_history == "aquatic" ,] # of those records, pull out only those that have aquatic life histories
  aq_count <- sum(list_aq$Abundance) # sum up all of the "presence" tallies for aquatic insects for that sample
  all_count <- sum(list$Abundance) # sum up total number of families found in sample
  per_aq <- (aq_count / all_count) # calculate the percent of each sample that is aquatic families
  occ_aq[i,1] <- sam # save the sample name
  occ_aq[i,2] <- per_aq # save the total aquatic relative abundance
}

# Name columns
names(occ_aq)[names(occ_aq) == "X1"] <- "sampleID"
names(occ_aq)[names(occ_aq) == "X2"] <- "percent_aquatic_occ"

# Add in sample information
aquatic <- merge(occ_aq, aquatic, by = "sampleID")

# Add in the number of reads per sample (sequencing depth) ---------------------

# Only keep the columns we need
depth_postprune <- depth_postprune[,c("sampleID", "TotalReads")]

# Merge with larger dataset
aquatic <- merge(aquatic, depth_postprune, by = "sampleID")

# Add column for "brood size at time of sampling" ------------------------------

aquatic$brood_size_sampling_time <- NA

for (i in 1:length(aquatic$sampleID)){
  if (aquatic$Age[i] == "6"){
    aquatic$brood_size_sampling_time[i] <- aquatic$Brood_Size_Day6[i]
  }
  if (aquatic$Age[i] == "12"){
    aquatic$brood_size_sampling_time[i] <- aquatic$Current_Brood_Size[i]
  }
  if (aquatic$Age[i] == "15"){
    aquatic$brood_size_sampling_time[i] <- aquatic$Current_Brood_Size[i]
  }
}

# Create a data frame which includes both nestling and adult female information together ---------
# Give each nestling its own line
# We will include adult females that don't have fecal samples here

# Pull out nestlings from aquatic 
aquatic_nestlings <- aquatic[aquatic$Adult_or_Nestling == "Nestling" ,] # Select nestlings

# Work in adult info from s_info file
captures <- s_info[s_info$Adult_or_Nestling == "Adult" ,]
captures <- captures[captures$Sex == "F" ,]

# Some third capture females do not have wing measurements, so import the wing measurements from their first capture
bands <- unique(captures$Individual_Band)

# The for loop that directly follows originally gave an error message because 
# there were some birds that were re-nests, but their first capture at their new
# box was listed as a first capture, and so they had multiple "first captures."
# In the following code, I rename those birds' "first" captures that actually 
# occurred at their re-nest nestbox as "renest"
# 278172145, cap 1 needing renest = 6/13/19
captures$Capture_Number[(captures$Individual_Band == "278172145" & captures$Capture_Date == "13-Jun-19")] <- "renest"
# 278172146, cap 1 needing renest = 6/9/19
captures$Capture_Number[(captures$Individual_Band == "278172146" & captures$Capture_Date == "09-Jun-19")] <- "renest"
# 278172231, cap 1 needing renest = 6/21/19
captures$Capture_Number[(captures$Individual_Band == "278172231" & captures$Capture_Date == "21-Jun-19")] <- "renest"
# 281128404, cap 1 needing renest = 6/17/19
captures$Capture_Number[(captures$Individual_Band == "281128404" & captures$Capture_Date == "17-Jun-19")] <- "renest"

# Make a for loop to put flatwing measurement from first capture as flatwing measurement for all captures
for (i in 1:length(bands)){
  band <- bands[i]
  flatwing <- captures$Flat_Wing[captures$Individual_Band == band & captures$Capture_Number == 1]
  cap2 <- captures[captures$Individual_Band == band & captures$Capture_Number == 2 ,]
  if (length(cap2$Individual_Band > 0)){
    captures$Flat_Wing[captures$Individual_Band == band & captures$Capture_Number == 2] <- flatwing
  }
  cap3 <- captures[captures$Individual_Band == band & captures$Capture_Number == 3 ,]
  if (length(cap3$Individual_Band > 0)){
    captures$Flat_Wing[captures$Individual_Band == band & captures$Capture_Number == 3] <- flatwing
  }
}

# Extract only third capture birds
captures <- captures[captures$Capture_Number == "3" ,] 

# Add prefix to all column names to delineate between adults and nestlings
colnames(captures) <- paste("ad", colnames(captures), sep = "_")
colnames(aquatic_nestlings) <- paste("n", colnames(aquatic_nestlings), sep = "_")

# Rename appropriate columns so they match and the data frames can be merged
captures <- dplyr::rename(captures,
                          Species = ad_Species,
                          Location = ad_Location,
                          Site = ad_Site,
                          Nest = ad_Nest,
                          Year = ad_Exp_Year,
                          site_box_year = ad_site_box_year,
)

aquatic_nestlings <- dplyr::rename(aquatic_nestlings,
                                   Species = n_Species,
                                   Location = n_Location,
                                   Site = n_Site,
                                   Nest = n_Nest,
                                   Year = n_Exp_Year,
                                   site_box_year = n_site_box_year,
)

aquatic_nestlingsadults <- merge(aquatic_nestlings, captures, all.x = FALSE)
# Note that some mothers were not captured a third time -- so their nestlings get
# taken out of this dataset.

# When it's available, we would like the adult samples info (i.e. percent aquatic)
# so that we can compare adult aquatic with nestling survival
aquatic_F_prov <- aquatic[aquatic$Sex == "F" & aquatic$Capture_Number == "3" ,]
aquatic_F_prov <- subset(aquatic_F_prov, select = c(Individual_Band, percent_aquatic_ra, percent_aquatic_occ, num_fam))
aquatic_F_prov <- dplyr::rename(aquatic_F_prov,
                                ad_Individual_Band = Individual_Band,
                                ad_percent_aquatic_ra = percent_aquatic_ra,
                                ad_percent_aquatic_occ = percent_aquatic_occ,
                                ad_num_fam = num_fam
)

# Add in this percent info
aquatic_nestlingsadults <- merge(aquatic_nestlingsadults, aquatic_F_prov, by = "ad_Individual_Band", all.x = TRUE)

# Write csvs to use in data analysis script ------------------------------------

write.csv(aquatic, here("2_modified_data/aquatic.csv"))
write.csv(aquatic_nestlingsadults, here("2_modified_data/aquatic_nestlingsadults.csv"))

# Percent aquatic with relative abundance vs. occurrence -- FIGURE S3 ----------
# Check on the relationship between percent aquatic when it is calculated
# with relative abundance vs. occurrence

p <- ggplot(aquatic, aes(x=percent_aquatic_ra, y = percent_aquatic_occ)) + geom_point() +
  geom_smooth(method=lm) + theme_classic() + xlab("Percent aquatic via relative abundance") +
  ylab("Percent aquatic via occurrence") +
  scale_y_continuous(limits = c(0,1)) +
  theme(axis.title=element_text(size = 18)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) # add a line showing a 1 to 1 relationship
ggsave(here("3_r_scripts/figs/figs_descriptive/Figure_S3_percent_aquatic_compare.png"), width = 5, height = 5, device = "png")

cor <- cor.test(aquatic$percent_aquatic_ra, aquatic$percent_aquatic_occ)
cor
# There is a low correlation between percent aquatic when calculated with relative
# abundance vs. occurrence.
