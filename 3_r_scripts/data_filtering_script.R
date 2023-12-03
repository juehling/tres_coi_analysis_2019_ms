# Data filtering for "Predictors and consequences of diet composition in a declining
# generalist aerial insectivore"

# Written by: Conor Taff and Jenny Uehling

# Last updated: 12/2/2023

# Run under R Studio with R version 4.3.1 on a Mac OS

# This code takes raw sequence data and metadata and converts it into formats
# useable in future analyses. It also performs basic calculations about the
# nature of the sequencing data and filters just to the data needed for this
# project.

# AMPtk processing pipeline ----------------------------------------------------

# Here, we describe the general pipeline for AMPtk that we used to process raw
# sequences. This section contains notes on how we did this in AMPtk, and the
# steps that brought us up to importing the files into R.

# All of the raw sequences for this project were processed using the `AMPtk`
# pipeline in the command line [@amptk]. We followed very closely the workflow
# described on the project website `amptk.readthedocs.io/en/latest`. Before
# running `AMPtk`, you need to install the software in a conda environment,
# install `USEARCH`, and download the COI database. Those steps are described in 
# detail on the website and not repeated here.

# A few notes about installation and setting up for the pipeline:
  
#  1. It seems that currently there is no way to run AMPtk on Mac with OS Catalina 
#     because it can't execute the 32bit `USEARCH` program. Use an older computer, 
#     virtual machine, or cloud computing.

#  2. The pipeline expects your sequences to be named: 
#     `sampleID_whatever_filler_R1/F1.fastq.gz`. The sequences from the Cornell
#     Biotechnology Resources Center (BRC) come back with the sample name in the 
#     middle. We haven't yet figured out how to make `AMPtk` recognize that, so 
#     we've been batch-renaming the files to get the sample name first.

# Once `AMPtk` is installed, the entire pipeline involves just three commands
# (will take several hours on a home laptop to run with samples in hundreds).
# Output files will be saved in the folder where you have set up your conda 
# environment. These commands are given below, along with descriptions of what
# they are doing.

# 1. Pre-processing sequence data: Takes a set prefix to add to all files and 
#    the forward and reverse primer sequences. Specify input folder where 
#    sequences are stored.

# > `amptk illumina -i /directory/with/sequences -o trescoi` 
# > `-f GCHCCHGAYATRGCHTTYCC -r TCDGGRTGNCCRAARAAYCA`

# 2. Denoising with `DADA2` [@dada2]: This takes in the processed sequences from
#    step 1 and applies the denoising algorithm that identifies ASVs and models 
#    sequencing error to arrive at a final list of ASVs for F/R reads. Forward 
#    and reverse reads are then merged together.

# > `amptk dada2 -i tres.coi.demux.fg.gz --platform illumina -o trescoi`

# 3. Assigning taxonomy: This uses the default suggestions from `AMPtk` 
#    documentation, but there are many other options. It applies the 'hybrid 
#    taxonomy algorithm', which looks at matches to i) Global Alignment to the 
#    downloaded COI arthropod database, ii) SINTAX classification, and iii) UTAX
#    classification. More information is on the website, but basically it retains
#    the best hit that also has the most levels of taxonomic ranks. The arthropod 
#    COI database needs to be downloaded and set up for this to run.

# > `amptk taxonomy -i trescoi.cluster.otu_table.txt -f trescoi.cluster.otus.fa -d COI`

# Those three commands produce a bunch of output files, but only three are needed
# to pull into `R` for subsequent analyses.

# - trescoi_2022_09r2.cluster.otu_table.txt has samples in columns and OTUs in rows with reads in each cell
# - trescoi_2022_09r2.cluster.taxonomy.txt has full taxonomic information for each OTU
# - trescoi_2022_09r2.mapping_file.txt has one row for each sample to add metadata

# All three of these files are read into the `R` script in this repository as the 
# basis of everything here with the addition of sample metadata.

# NOTE: There are many other optional steps in `AMPtk`. For example, one set
# of commands works to try to remove barcode bleed over where sequences from one
# barcode/sample are accidentally attributed to other samples. There are 
# different types of classifiers and options, etc. What we've done here is
# the basic or bare-bones pipeline.

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

# Load and wrangle tree swallow data -------------------------------------------

# The tree swallow metadata come directly from the Vitousek Lab Tree Swallow 
# Database, which is maintained by the Vitousek Lab Manager in Microsoft Access.

# Downloaded Excel file with individual tree swallow metadata 
# (Captures_Hormone_Bleeding_Blood_DNA_12.09.21.xlsx) and converted it to .csv

# This file is housed outside of the Github because it contains all records
# from the Vitousek lab long-term dataset; in the following code, we extract
# just the data we need and save it in this Github repository.
s_info <- read.csv("~/Dropbox/tres_database_data/Captures_Hormone_Bleeding_Blood_DNA_12.09.21.csv")

# Delete all rows not from 2019.
s_info <- s_info[!is.na(s_info$Exp_Year), ] # First must delete single row with NA for year.
s_info <- s_info[s_info$Exp_Year == "2019" ,]

# Delete columns that are empty or not needed to prevent dataframe from being too unwieldy.
s_info <- subset(s_info, select = -c(Tail, Fat, Culmen, Ectoparasites_Presence, Genetic_Mom, 
                                     Genetic_Dad, Natal_Nest, Trap_set, Disturbance_Time, Bleed1_Latency,
                                     Bleed1_Amount, Bleed1_Plasma_nb, Base_Plasma_Sample_Box,
                                     Baseline_Sample_Location, Bleed1_CORT_Raw, Bleed1_OXY, Bleed1_dROMs,
                                     BUTY, Bleed2_Time, Bleed2_Amount, Bleed2_Plasma_nb, Stress_Plasma_Sample_Box,
                                     Stress_Sample_Location, Bleed2_CORT_Raw, Bleed2_OXY, Bleed2_dROMs,
                                     Inject_1_Time, Inject_2_Time, Bleed3_Time, Bleed3_Amount,
                                     Bleed3_Plasma_nb, Dex_Plasma_Sample_Box, Dex_Sample_Location,
                                     Bleed3_CORT_Raw, Release_Time, Feathers_Sample_nb, Telomeres,
                                     Telomere_Sample_nb, Telomeres_Box, Telomeres_Sample_Location, Microbiome,
                                     Microbiome_Sample_nb, Microbiome_Box, Microbiome_Sample_Location, Bolus_nb,
                                     Blood_Smear_Collected, Blood_DNA_Collected, Bacteria_Killing_Assay,
                                     OXY_AssayCode, dROMS_AssayCode, Manipulated_after_1st_capture,
                                     Method, Leg_Banded, Status, Measured_by, Blood_DNA_Sample_nb, Lifetag_nb,
                                     RIA_or_EIA, CORT_AssayCode, Extraction_Efficiency, Bleed1_Glucose, Bleed2_Glucose,
                                     Bleed3_Glucose, RFID_On, Individual_RFID, Stress_Series_Type, Bleed1_CORT_corrected_value,
                                     Bleed2_CORT_corrected_value, Inject_1_Type, Inject_1_Dose, Inject_2_Type, Inject_2_Dose, 
                                     Bleed3_CORT_corrected_value))

# Create column for unit_box_year to use to match up treatments with nestlings/nests
s_info$site_box_year <- paste(s_info$Site, s_info$Nest, s_info$Exp_Year, sep="_")

# Import file with information not contained in individual tree swallow metadata:
# day 6 nestling sample metadata and blanks metadata (i.e., "blanks" from DNA extraction)
extra_info <- read.csv(here("1_raw_data", "extra_samples.csv"))

# Create a column in extra_info that has the same information as site_box_year
extra_info$site_box_year <- paste(extra_info$Site, extra_info$Nest, extra_info$Exp_Year, sep="_")

# Delete "cap_doy" column from extra_info because it is not needed
extra_info <- subset(extra_info, select = -c(cap_doy))

# Reclassify columns in s_info and extra_info so they bind
s_info$Individual_Band <- as.character(s_info$Individual_Band)
s_info$Exp_Year <- as.character(s_info$Exp_Year)
s_info$Capture_Number <- as.character(s_info$Capture_Number)
s_info$Age <- as.character(s_info$Age)
extra_info$Age <- as.character(extra_info$Age)

# Combine s_info and extra_info
s_info <- bind_rows(s_info, extra_info)

# Downloaded Excel file with tree swallow nest metadata 
# (Nest_Records 12.09.21.xlsx) and converted it to .csv

# This file is housed outside of the Github because it contains all records
# from the Vitousek lab long-term dataset; in the following code, we extract
# just the data we need and save it in this Github repository.
nest_info <- read.csv("~/Dropbox/tres_database_data/Nest_Records_12.09.21.csv")

# Delete all rows not from 2019.
nest_info <- nest_info[nest_info$Exp_Year == "2019" ,]

# Create a column in nest_info that has the same information as site_box_year
nest_info$site_box_year <- paste(nest_info$Site, nest_info$Nest, nest_info$Exp_Year, sep="_")

# Create a column in nest_info that has the average nestling mass on day 6
nest_info$avg_day6_mass <- nest_info$Day6_Brood_Mass/nest_info$Brood_Size_Day6

# Delete the column that's in "s_info" that has experiment and treatment.
# This is because we're going to use all the experiment info from the nest spreadsheet.
s_info <- subset(s_info, select=-c(Experiment, Individual_Treatment))

# Now, we'll extract just the columns we need from nest_info
nest_info <- subset(nest_info, select=c(Nest_Experiment, Nest_Treatment, 
                                        site_box_year, Day6_Brood_Mass, Brood_Size_Day6, Nest_Fate, avg_day6_mass))

# Put nest information in "s_info"
s_info <- left_join(s_info, nest_info, by = "site_box_year")

# Create a column that lists both the experiment AND the treatment
s_info$exp_treat <- paste(s_info$Nest_Experiment, s_info$Nest_Treatment, sep="_")

# Correctly classify control birds -- there are different naming conventions across individuals
table(s_info$exp_treat) # Figure out all of the "control" treatments listed
s_info$exp_treat[s_info$exp_treat == "_"] <- "Control"
s_info$exp_treat[s_info$exp_treat == "_Control"] <- "Control"
s_info$exp_treat[s_info$exp_treat == "NA_NA"] <- "Control"
table(s_info$exp_treat)

# Rename fecal sample number to sampleID
s_info <- s_info %>% dplyr::rename(
  sampleID = Fecal_Sample_nb
)

# Save "s_info" for future use -- it will be our metadata file.
write.csv(s_info, "2_modified_data/s_info.csv")

# Load and wrangle sequencing data ---------------------------------------------

# The prefix used for AMPtk processing is tres_coi.
# This is set in AMPtk and all files produced there have this prefix.
amptk_prefix <- "trescoi"

# Load the number of reads by taxa per sample table. Format for phyloseq.
otu_ab <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2022_09r2.cluster.otu_table.txt")))

# P19N379 is mislabeled; it should be labelled as P19N479.
otu_ab <- dplyr::rename(otu_ab, F05xP19N479 = F05xP19N379)

rownames(otu_ab) <- otu_ab$X.OTU.ID   # give row names from otu names
otu_ab <- otu_ab[, 2:ncol(otu_ab)]    # remove the column of otu names
  
# Read the mapping table
# This 'mapping' table from AMPtk is mostly blank but has all the sample
# names so it can be joined to actual sample metadata. It's also possible
# to merge in the sample metadata within the AMPtk pipeline.
map <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2022_09r2.mapping_file.txt")))
  
# P19N379 is mislabeled; it should be labelled as P19N479.
map$X.SampleID <- gsub("F05xP19N379", "F05xP19N479", map$X.SampleID)
map$phinchID <- gsub("F05xP19N379", "F05xP19N479", map$phinchID)

# Make sampleID its own column
for(i in 1:nrow(map)){
  map$sampleID[i] <- strsplit(map$X.SampleID[i], "x")[[1]][2]
}
  
# write.table(map, "map.txt", sep = "\t") # to save a copy of the mapping file
  
# Read the otu taxonomy table
otu_tax <- read.delim(here("1_raw_data", paste0(amptk_prefix, "_2022_09r2.cluster.taxonomy.txt")))
rownames(otu_tax) <- otu_tax$X.OTUID

# The taxonomy result from AMPtk is in one long string of text. This is splitting up the string
# and filling in a bunch of different columns. 
for(i in 1:nrow(otu_tax)){
  temp <- otu_tax$taxonomy[i]
  temp2 <- strsplit(temp, "\\|")[[1]][3]
  temp3 <- strsplit(temp2, ":")[[1]][2]
  otu_tax$search_hit[i] <- strsplit(temp, "\\|")[[1]][1]
  otu_tax$hit_score[i] <- strsplit(temp, "\\|")[[1]][2]
  otu_tax$database[i] <- strsplit(temp2, ":")[[1]][1]
  otu_tax$accession[i] <- strsplit(temp3, ";")[[1]][1]
  otu_tax$kingdom[i] <- strsplit(strsplit(temp2, "k:")[[1]][2], ",")[[1]][1]
  otu_tax$phylum[i] <- strsplit(strsplit(temp2, "p:")[[1]][2], ",")[[1]][1]
  otu_tax$class[i] <- strsplit(strsplit(temp2, "c:")[[1]][2], ",")[[1]][1]
  otu_tax$order[i] <- strsplit(strsplit(temp2, "o:")[[1]][2], ",")[[1]][1]
  otu_tax$family[i] <- strsplit(strsplit(temp2, "f:")[[1]][2], ",")[[1]][1]
  otu_tax$genus[i] <- strsplit(strsplit(temp2, "g:")[[1]][2], ",")[[1]][1]
  otu_tax$species[i] <- strsplit(strsplit(temp2, "s:")[[1]][2], ",")[[1]][1]		
  }
  
# Replace database mismatches caused by matches that aren't from BOLD records
otu_tax$database <- gsub("None;k", "None", otu_tax$database)
otu_tax$accession <- gsub("Animalia,p", "None", otu_tax$accession)
  
# For phyloseq this has to be added as a matrix rather than data frame    
otu_tax <- as.matrix(otu_tax)

# This is saving just the taxonomic ranks rather than the database info.
otu_tax_2 <- otu_tax[, 10:16]
  
# Add in the sample information to each sample
map_2 <- join(map, s_info, "sampleID", "left", "first")
rownames(map_2) <- map_2$X.SampleID
  
# Print out some summary information about this dataset
# This will identify how many samples there are in each category of site,
# age, and capture number.
summary_table <- map_2 %>% dplyr::count(Site, Age, Capture_Number)

# This will identify samples that don't match the metadata file and write them as a separate output
missing <- subset(map_2, is.na(map_2$Individual_Band) == TRUE) # none are missing
# write.table(missing, "missing_info.txt", sep = "\t")

# Create metadata file to upload to NCBI ---------------------------------------

# Save list of samples
sample_list <- as.data.frame(map_2$sampleID)

# Rename column so it can merge with s_info
sample_list <- dplyr::rename(sample_list, "sampleID" = "map_2$sampleID")

# Merge with s_info
sample_list <- left_join(sample_list, s_info, by = "sampleID")

# Rename columns to match with NCBI format requested
sample_list <- dplyr::rename(sample_list, "sample_name" = "sampleID")

# Make column for host_subject_ID for NCBI upload -- requires a column
# independent of sampleID that is unique to each sample
sample_list$NCBI_host_subject_ID <- NA

for (i in 1:length(sample_list$sample_name)){
  if (sample_list$Adult_or_Nestling[i] == "Nestling"){
    if (sample_list$Age[i] == "6"){
      sample_list$NCBI_host_subject_ID[i] <- paste("Day_6_nestling", sample_list$site_box_year[i], sample_list$sample_name[i], sep = "_")
    }
    if (sample_list$Age[i] == "12"){
      sample_list$NCBI_host_subject_ID[i] <- paste("Day_12_nestling", sample_list$Individual_Band[i], sample_list$Age[i], sep = "_")
    }
    if (sample_list$Age[i] == "15"){
      sample_list$NCBI_host_subject_ID[i] <- paste("Day_15_nestling", sample_list$Individual_Band[i], sample_list$Age[i], sep = "_")
    }
  }
  if (sample_list$Adult_or_Nestling[i] == "Adult"){
    sample_list$NCBI_host_subject_ID[i] <- paste("Adult", sample_list$Individual_Band[i], sample_list$Capture_Date[i], sep = "_")
  }
  if (sample_list$Adult_or_Nestling[i] == "neg_control"){
    sample_list$NCBI_host_subject_ID[i] <- paste("Blank", sample_list$sample_name[i], sep = "_")
  }
}

# Save file for export to NCBI metadata file
write.csv(sample_list, "2_modified_data/metadata_for_NCBI.csv")

# Build initial phyloseq object ------------------------------------------------

TAX = tax_table(otu_tax_2)
SAM = sample_data(map_2)
OTU = otu_table(otu_ab, taxa_are_rows = TRUE)

coi_ps <- phyloseq(OTU, TAX, SAM)

# Save initial phyloseq object
saveRDS(coi_ps, "2_modified_data/coi_ps.rds")

# Filter out other samples not needed ------------------------------------------

# Take out samples that we shouldn't use for further analyses due to labeling issue
coi_ps2 <- subset_samples(coi_ps, sampleID != "P19N009")
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N379")
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N478")
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N479") # This sample's labeling was corrected earlier, but it's still unclear if the incorrect sample was extracted.
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N496")
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N497")
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N498")
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N499")

# Take out adult male samples from Unit 4 and Turkey Hill because there are not
# enough samples to conclude anything, so we can't use them in this project.
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N215") # M 4.80
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N268") # M 4.93
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N366") # M 4.61
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N184") # M TH.3
coi_ps2 <- subset_samples(coi_ps2, sampleID != "P19N257") # M TH.86

# Summarize sequences before filtering -----------------------------------------

# Blanks: 26

# October 2019 sequencing run: 100 (NOT including blanks)
  # Plate 3: 5
  # Plate 4: 91 / 4 blanks
  # Plate 5: 4
# First July 2020 sequencing run: 335 (NOT including blanks)
  # Plate 1: 92 / 4 blanks
  # Plate 2: 92 / 4 blanks
  # Plate 3: 91 / 5 blanks
  # Plate 4: 60 / 3 blanks
# Second July 2020 sequencing run: 21
  # Plate 1: 21
# November 2020 sequencing run had: 105 (NOT including blanks or captive birds)
  # Plate 1: 92 / 4 blanks
  # Plate 2: 13 / 2 blanks
# Total: 561
# Total plus blanks: 587
# Deleted samples due to issues (explained above): 13
# Total samples that should be in dataset: 574 (excluding blanks: 548). Total that are: 573 (excluding blanks: 547)
# P19N035 is missing from the sequencing results; it was sequenced in July (second July 2020 run, on Plate 1).
# It had one DemuxRead, but it was removed sometime during the AMPtk processing.
# Notes from extraction said there was very little initial sample. Determined to be a sequencing fail.

# Create a record of the sequencing depth of each sample before filtering to arthropods
depth_no_arth_filter <- data.frame(as(sample_data(coi_ps2), "data.frame"),
                                   TotalReads = sample_sums(coi_ps2), keep.rownames = TRUE)

  ## Calculate number of adults and nests sampled ------------------------------

  # Determine number of adults sampled
  length(unique(depth_no_arth_filter$Individual_Band[depth_no_arth_filter$Adult_or_Nestling == "Adult"]))

  # Determine number of nests sampled
  nestlings <- depth_no_arth_filter[depth_no_arth_filter$Adult_or_Nestling == "Nestling" ,]
  length(unique(nestlings$site_box_year))

  ## Negative controls ---------------------------------------------------------

  # Extract negative controls for summary statistics
  depth_neg_controls <- depth_no_arth_filter[depth_no_arth_filter$Adult_or_Nestling == "neg_control" ,]
  
  # Sum up total number of reads
  sum(depth_neg_controls$TotalReads)

  # Determine the range of reads across samples
  range(depth_neg_controls$TotalReads)

  # Determine mean and median number of reads per sample
  mean(depth_neg_controls$TotalReads)
  median(depth_neg_controls$TotalReads)
  
  # Examine ASVs in blanks
  coi_ps2_blanks <- subset_samples(coi_ps2, Adult_or_Nestling == "neg_control")
  (coi_ps2_blanks = prune_taxa(taxa_sums(coi_ps2_blanks) > 0, coi_ps2_blanks)) # take out taxa from tax_table that aren't found in blanks
  tax_tab_examine_blanks <- data.frame(tax_table(coi_ps2_blanks))
  
  # Determine number of ASVs
  length(tax_tab_examine_blanks$kingdom)

  ## Fecal samples -------------------------------------------------------------

  # Take out negative controls for summary statistics
  depth_no_arth_filter_negcontrolremoved <- depth_no_arth_filter[depth_no_arth_filter$Exp_Year != "neg_control" ,]

  # Sum up total number of reads
  sum(depth_no_arth_filter_negcontrolremoved$TotalReads)

  # Determine the range of reads across samples
  range(depth_no_arth_filter_negcontrolremoved$TotalReads)

  # Determine mean and median number of reads per sample
  mean(depth_no_arth_filter_negcontrolremoved$TotalReads)
  median(depth_no_arth_filter_negcontrolremoved$TotalReads)

  depth_no_arth_filter <- subset(depth_no_arth_filter_negcontrolremoved, select = c(sampleID, Individual_Band, Age, Site, Nest, TotalReads)) # for readability
  write.csv(depth_no_arth_filter, "2_modified_data/depth_no_arth_filter.csv")

  # Examine ASVs in samples (exclude blanks)
  coi_ps2_noblanks <- subset_samples(coi_ps2, Exp_Year != "neg_control")
  (coi_ps2_noblanks = prune_taxa(taxa_sums(coi_ps2_noblanks) > 0, coi_ps2_noblanks)) # take out taxa from tax_table that aren't found in samples (i.e., that are only found in blanks)
  tax_tab_examine <- data.frame(tax_table(coi_ps2_noblanks))

  # Determine number of ASVs
  length(tax_tab_examine$kingdom)

# Filter to arthropods ---------------------------------------------------------

# Subset just to arthropods
coi_ps2 <- subset_taxa(coi_ps2, phylum == "Arthropoda")

# Extract the unique arthropod families found in the dataset
unique_families <- get_taxa_unique(coi_ps2, taxonomic.rank = "family")

# Save file with list of unique families to use to research aquatic vs.
# terrestrial families.
write.csv(unique_families, here("4_other_outputs/unique_families.csv"))
# This file was then taken out of R to research aquatic vs. terrestrial
# families, and will be re-imported later.

# Check the number of reads across fecal samples and blanks -- FIGURE S2 -------

# Histogram of sample reads. This is counting a sum of arthropod reads for each sample.
depth <- data.frame(as(sample_data(coi_ps2), "data.frame"),
                    TotalReads = sample_sums(coi_ps2), keep.rownames = TRUE)

depth$Adult_or_Nestling[depth$Adult_or_Nestling == "neg_control"] <- "Negative control"

p <- ggplot(depth, aes(log(TotalReads))) + geom_histogram(fill = "slateblue") + 
  ylab("Number of samples") + xlab("log(Reads)") +
  theme_classic() + geom_vline(xintercept = log(100), linetype = "dashed", col = "coral3", size = 1) + 
  geom_text(x = log(100) - 0.5, y = 30, label = "100 Reads", angle = 90)

p2 <- p + facet_grid(~ Adult_or_Nestling) +    # same but splitting out adult/nestling/negative_control
      theme(strip.text.x = element_text(size = 12))

# Save the histograms
ggsave(here("3_r_scripts/figs/figs_descriptive/total_reads.png"), plot = p, width = 8, height = 4.5, device = "png")
ggsave(here("3_r_scripts/figs/figs_descriptive/Figure_S2_total_reads_split.png"), plot = p2, width = 8.2, height = 4, device = "png")

# Extract information about blanks ---------------------------------------------

negative_controls <- depth[depth$Adult_or_Nestling == "Negative control" ,]

# Sum up total number of reads after filtering to arthropods
sum(negative_controls$TotalReads)

# Determine the range of reads across samples
range(negative_controls$TotalReads)

# Determine mean and median number of reads per sample
mean(negative_controls$TotalReads)
median(negative_controls$TotalReads)

# Remove blanks (negative controls) --------------------------------------------

# Remove negative controls
coi_ps2 <- subset_samples(coi_ps2, Age != "neg_control")

# Save phyloseq object ---------------------------------------------------------

saveRDS(coi_ps2, "2_modified_data/coi_ps2.rds")
