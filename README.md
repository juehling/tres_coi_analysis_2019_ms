
# Purpose

This R project contains all data and R scripts used to produce the analyses, figures, and tables in "Predictors and consequences of diet composition in a declining generalist aerial insectivore." The project takes sequences from tree swallow fecal samples, identifies their contents, and explores predictors and consequences of diet composition. What follows are descriptions of each folder and how to use them to reproduce the analyses for this project.

Authors: Jennifer J. Uehling, Conor C. Taff, Jennifer L. Houtz, Paige M. Becker, Allison S. Injaian, and Maren N. Vitousek

# Folders

## 1_raw_data

Contains raw data for this project.

- **Files with the prefix "trescoi_2022_09r2"**: results from running the sequencing result fastq files through the AMPtk pipeline, including the otu table, taxonomy file, and mapping file. See *data_filtering_script.R* for more information about the AMPtk pipeline.
- **extra_samples.csv**: contains records of samples sequenced that do not have metadata in the Vitousek Lab database files. This includes blanks and day 6 nestlings that do not yet have individual identification (they were not marked with uniquely-number bands until day 12).

The raw sequencing files (fastq files) for this project, which were run through the AMPtk pipeline, are uploaded to NCBI’s Sequence Read Archive (Bioproject: PRJNA884756).

## 2_modified_data

Contains files created with R scripts (contained in 3_r_scripts). All files in this folder are created in the r scripts listed below. The metadata file for this project is *s_info.csv*, and it is created in data_filtering_script.R. The metadata file comes from the larger Vitousek Lab database files, which are stored locally on the J. Uehling's computer. The meaning of each column in *s_info.csv* is listed at the end of this page.

## 3_r_scripts

Contains all scripts used for data filtering and analysis. The scripts should be run in the order listed. Outline of scripts are as follows:

- **data_filtering_script.R**: Takes files created with AMPtk pipeline (explained in script), appends metadata (*s_info.csv*), filters to arthropod, explores families, checks sample effort, removes negative controls, and saves all data in phyloseq object.
- **relative_abundance_occurrence_calculations**: Performs more filtering steps. Calculates general summary statistics about fecal sample contents and creates summary figures. Calculates the relative abundance of aquatic insects in each fecal sample, using relative abundance of aquatic insect distinct taxonomic group reads, occurrence of aquatic insect distinct taxonomic group reads, and occurrence of aquatic insect family reads.
- **diversity_calculations.R**: Performs filtering steps. Calculates diversity indices for each sample.
- **data_analysis_script.R**: Performs modeling and makes tables and figures with products created in three aforementioned R scripts.

**Folder "model_outputs"**: Tables from models, produced in *data_analysis_script.R*

**Folder "figs"**: Figures. Includes figures created in *relative_abundance_occurrence_calculations.R* and *data_analysis_script.R*

## 4_other_outputs

Contains spreadsheets with the unique arthropod families identified in the fecal samples plus information about their immature stages (i.e., aquatic, terrestrial, both, unknown).

- *unique_families.csv*: created in *data_filtering_script.R*. Contains a list of the unique families found in the fecal samples.
- *unique_families_aquatic_terrestrial_IDs.csv*: populated outside of R. Contains identifications of the habitat of the larval stage of each insect family (aquatic, terrestrial, both, or unknown). Columns are as follows:
	- *order*: order
	- *family*: family
	- *life_history*: the habitat for the larval stage of the family
	- *common_name*: common name of arthropod family
	- *sources*: sources of information about the arthropod family
	- *Notes*: notes about arthropod family's larval stage.

# Column names of raw data files (stored in *1_raw_data*) and metadata file (*s_info.csv*)

Below are the column names for each of the raw data files include in the analysis.

## extra_samples.csv

- *Fecal_Sample_nb*: fecal sample ID number
- *Species*: species of bird, using the USGS four letter banding code (tree swallow, or TRES)
- *Location*: city location of capture
- *Site*: site of capture
- *Nest*: nest box number
- *Exp_Year*: year of capture
- *Adult_or_Nestling*: identifies whether the bird was an adult or a nestling
- *Capture_Date*: date the bird was captured (day-month-year, in format dd-Mon-yy)
- *cap_doy*: the Julian day of the year when the bird was captured
- *Age*: age of bird in days

## trescoi_2022_09r2.mapping_file.txt

- *SampleID*: sample well location in plate from when samples were submitted to the Cornell BRC and fecal sample ID number, with an "x" in between
- *BarcodeSequence*: barcode assigned to each sample by Cornell BRC
- *LinkerPrimerSequence*: forward primer (BF2) sequence
- *RevBarcodeSequence*: not applicable, residual column from AMPtk pipeline
- *ReversePrimer*: reverse primer (BR2) sequence
- *phinchID*: duplicate of *SampleID* created by AMPtk pipeline, not applicable
- *DemuxReads*: number of reads in sample after demultiplexing and before denoising
- *Treatment*: not applicable, residual column from AMPtk pipeline

## trescoi_2022_09r2.cluster.otu_table.txt

- Each column is a sampleID and each row is an ASV. The numbers in each cell indicate the number of reads of each ASV (row) in each sample (column).

## trescoi_2022_09r2.cluster.taxonomy.txt

- *OTUID*: ASV ID number (assigned by AMPtk)
- *taxonomy*: The taxonomic identification of the ASV. AMPtk outputs this as a long string of text and includes some additional information in the output; see AMPtk documentation for more information.
- *USEARCH*: Taxonomic identification of the ASV from the downloaded COI arthropod database (see notes in *data_filtering_script.R*)
- *SINTAX*: Taxonomic identification of the ASV from SINTAX
- *UTAX*: Taxonomic identification of the ASV from UTAX

## s_info.csv

- *Individual_Band*: USGS-issued band number of bird
- *Exp_Year*: year of capture
- *Adult_or_Nestling*: identifies whether the bird was an adult or a nestling
- *Sex*: sex of bird (male, female, or unknown)
- *Capture_Date*: date the bird was captured (day-month-year, in format dd-Mon-yy)
- *Species*: species of bird, using the USGS four letter banding code (tree swallow, or TRES)
- *Location*: city location of capture
- *Site*: site of capture
- *Nest*: nest box number
- *Capture_Number*: capture number of bird that year (1, 2, or 3)	
- *Age*: age of bird. Adult age possibilities are SY (second year), ASY (after second year), or AHY (after hatch year). Nestling age possibilities are 6, 12, and 15, in days.
- *Mass*: mass in grams
- *Bill.Head*: head plus bill length in mm
- *Flat_Wing*: flat wing length in mm
- *Current_Brood_Size*: brood size at time of capture
- *Nestling_Fate*: eventual fate of nestling (died, fledged, predated, unknown)
- *sampleID*: fecal sample ID number
- *Notes*: notes about the capture
- *site_box_year*: the site of breeding (Turkey_Hill, Unit_1, Unit_2, or Unit_4), the nest box, and the year, all separated by underscores
- *Nest_Experiment*: experiment that was occurring at the nest of the individual captured
- *Nest_Treatment*: specific treatment that was occurring at the nest of the individual captured
- *Day6_Brood_Mass*: mass of all nestlings in that nest on day 6 of nestling development
- *Brood_Size_Day6*: brood size on day 6 of nestling development
- *Nest_Fate*: eventual fate of nest (died, fledged, predation, abandoned)
- *avg_day6_mass*: average mass of each nestling in the nest, calculated by dividing Day6_Brood_Mass by Brood_Size_Day6
- *exp_treat*: Nest_Experiment and Nest_Treatment separated by an underscore
