library(tidyverse)
library(stringr)
library(Gviz)
library(ggpubr)
library(regioneR)

# Set WD ----------------------------------------------------------------------------------------------

# OMIM_wd <- Sys.getenv("OMIM_wd")
# setwd(OMIM_wd)

#setwd("/home/dzhang/projects/OMIM_wd")

# Load data -------------------------------------------------------------------------------------------

# Functions -------------------------------------------------------------------------------------------

source("~/R/Proof_Concept/visualise_ER_example_gviz_sonia.R")

# Main ------------------------------------------------------------------------------------------------

visualise_ER_example(ERs_w_annotation_df = head(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific,1000), 
                     txdb = ensembl_grch38_v92_genes_txdb, 
                     ensembl_gene_id = "ENSG00000145335",
                     tissues_to_plot = "frontalcortexba9", 
                     genome_build = "hg38",
                     gtex_split_read_table_mean_cov_df,
                     tissue_optimal_cut_off_max_gap_df,
                     get_constraint = F,
                     get_conserv = F,
                     propor_samples_split_read = 0.05,
                     extend_region_to_plot = 10000,
                     collapseTranscripts = "meta",
                     transcriptAnnotation = "gene",
                     aceview_annot = NULL,
                     add_custom_annot_track = NULL,
                     all_split_reads = F)

dev.print(file = "~/R/Proof_Concept/SNCA_OMIM_reannot_example.png", device = png, res = 600, width = 10, height = 11.69/2, units = "in")

visualise_ER_example(ERs_w_annotation_df = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, 
                     txdb = ensembl_grch38_v92_genes_txdb, 
                     ensembl_gene_id = "ENSG00000134086",
                     tissues_to_plot = c("brain_cerebellar_hemisphere", "lung"), 
                     genome_build = "hg38",
                     gtex_split_read_table_mean_cov_df,
                     tissue_optimal_cut_off_max_gap_df,
                     get_constraint = F,
                     get_conserv = F,
                     propor_samples_split_read = 0.05,
                     extend_region_to_plot = 3000,
                     collapseTranscripts = "meta",
                     transcriptAnnotation = "gene",
                     aceview_annot = NULL,
                     add_custom_annot_track = NULL,
                     all_split_reads = F)

dev.print(file = "~/R/Proof_Concept/ERLIN1_OMIM_reannot_example.png", device = png, res = 600, width = 10, height = 11.69/2, units = "in")
