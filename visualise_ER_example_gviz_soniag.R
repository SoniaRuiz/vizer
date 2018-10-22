library(tidyverse)
library(stringr)
library(Gviz) 
library(ggpubr)
library(regioneR)

# Set WD ----------------------------------------------------------------------------------------------

OMIM_wd <- "/home/dzhang/projects/OMIM_wd"
setwd(OMIM_wd)

# Load data -------------------------------------------------------------------------------------------

load(file = "results/annotate_ERs/ERs_optimised_cut_off_max_gap_all_tissues_w_annot_df.rda")

tissue_optimal_cut_off_max_gap_df <- read_delim("results/optimise_derfinder_cut_off/exon_delta_details_optimised_maxgap_cutoff.csv", delim = ",")

gtex_tissue_name_formatting <- read_delim("raw_data/gtex_tissue_name_formatting/OMIM_gtex_tissue_name_formatting.csv", delim = ",")

ensembl_gene_id_to_symbol_df_v92 <- read_delim("raw_data/data_description_files/ensembl_gene_id_to_symbol_df_v92.csv", delim = ",")

# Functions -------------------------------------------------------------------------------------------

##### First level #####

source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/generate_txDb_from_gtf.R")
source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/constraint/constraint_general/constraint_general_functions.R")

get_gtex_split_read_table_mean_cov_n_samples_df <- function(gtex_tissue_name_formatting){
  
  gtex_split_read_table_annotated_paths <- 
    list.files("/data/recount/GTEx_SRP012682/gtex_split_read_table_annotated_rda/", full.names = T)
  
  gtex_split_read_table_df <- 
    data_frame(gtex_split_read_table_annotated_paths = gtex_split_read_table_annotated_paths, 
               tissue = 
                 gtex_split_read_table_annotated_paths %>% 
                 str_replace("/.*/", "") %>% 
                 str_replace("_split_read_table_annotated.rda", "") %>% 
                 str_replace("brain", "") %>% 
                 str_replace_all("_", ""))
  
  gtex_mean_cov_df <- 
    data_frame(gtex_mean_cov_paths = 
                 list.files("/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/", full.names = T), 
               tissue = 
                 list.files("/data/recount/GTEx_SRP012682/gtex_mean_coverage/by_tissue_smfrze_use_me/", full.names = F) %>% 
                 str_replace("brain", "") %>% 
                 str_replace_all("_", ""))

  gtex_split_read_table_mean_cov_df <- 
    gtex_split_read_table_df %>% 
    inner_join(gtex_mean_cov_df) %>% 
    left_join(gtex_tissue_name_formatting, by = c("tissue" = "gtex_tissue_name_simplified"))
  
  return(gtex_split_read_table_mean_cov_df)
  
}

#' @param ERs_w_annotation_df dataframe with ER coordinates and annotation 
#' @param txdb txdb object with the definitions of genes inside 
#' @param ensembl_gene_id_to_symbol_df df with all the ensembl gene ids matched to symbols using biomart
#' @param gene_id gene id to plot 
#' @param tissues_to_plot tissues to plot
#' @param gtex_split_read_table_mean_cov_df df that maps the tissue to it's split read or mean coverage path
get_ER_table_to_display <- function(ERs_w_annotation_df, txdb, ensembl_gene_id_to_symbol_df, gene_id, 
                                    tissues_to_plot, gtex_split_read_table_mean_cov_df, extend_region_to_plot){
  
  ##### convert symbol to ENSG ID #####
  
  if(str_detect(gene_id, "ENSG")){
    
    ensembl_gene_id <- gene_id
    
  }else{
    
    ensembl_gene_id_to_symbol_df_dup <- 
      ensembl_gene_id_to_symbol_df %>%
      filter(!duplicated(external_gene_name), # removing instances whereby 2 ensembl ids match onto 1 
             external_gene_name == gene_id)
    
    if(nrow(ensembl_gene_id_to_symbol_df_dup) == 0){
      
      stop(str_c(gene_id, " not found in the ensembl database (v92)"))
      
    }else{
      
      ensembl_gene_id <- ensembl_gene_id_to_symbol_df_dup$ensembl_gene_id
      
    }
    
  }
  
  ensembl_grch38_v92_genes_txdb_genes <- genes(ensembl_grch38_v92_genes_txdb)
  
  gene_cord <- ensembl_grch38_v92_genes_txdb_genes[ensembl_grch38_v92_genes_txdb_genes$gene_id == ensembl_gene_id]
  seqnames_to_plot <- as.character(seqnames(gene_cord))
  start_to_plot <- (start(gene_cord) - extend_region_to_plot)
  end_to_plot <- (end(gene_cord) + extend_region_to_plot)
  bp_plotted <- end_to_plot - start_to_plot
  
  ERs_w_annotation_df_to_display <- 
    ERs_w_annotation_df %>% 
    filter(seqnames == seqnames_to_plot, start >= start_to_plot, end <= end_to_plot, 
           tissue %in% tissues_to_plot) %>% 
    dplyr::select(ER_seqnames = seqnames, ER_start = start, ER_end = end, ER_width = width, tissue, mean_coverage = value, ensembl_grch38_v92_region_annot, overlap_any_gene_v92_name, 
                  split_read_annotation_type = annotationType_split_read_annot, annotationType_split_read_annot, split_read_to_any_gene = uniq_genes_split_read_annot, 
                  mean_CDTS_percentile, mean_phast_cons_7, mean_phast_cons_100)
  
  return(ERs_w_annotation_df_to_display)
  
}

#' @param ERs_w_annotation_df dataframe with ER coordinates and annotation 
#' @param txdb txdb object with the definitions of genes inside 
#' @param ensembl_gene_id_to_symbol_df df with all the ensembl gene ids matched to symbols using biomart
#' @param gene_id gene id to plot 
#' @param tissues_to_plot tissues to plot
#' @param genome_build genome build - defaults to hg38 
#' @param gtex_split_read_table_mean_cov_df df that maps the tissue to it's split read or mean coverage path
#' @param get_constraint lgl vector indicating whether to plot constraint scores
#' @param get_conserv lgl vector indicating whether to plot conserv scores 
#' @param propor_samples_split_read percentage of samples required that have the split read to plot
#' @param extend_region_to_plot number of bps either side of the gene to add to around plot
#' @param aceview_annot defaults to NULL otherwise, takes the aceview txdb 
#' @param add_custom_annot_track here gives the option to input your own positions to plot as a GRanges
#' @param all_split_reads lgl vector determining whether to plot all split reads in the region, not just those that connect to ERs
#' @return Gviz plot to visualise the ERs in the genomic region or by gene
visualise_ER_example <- function(ERs_w_annotation_df, txdb, ensembl_gene_id_to_symbol_df, gene_id, tissues_to_plot, 
                                 genome_build = "hg38", gtex_split_read_table_mean_cov_df, 
                                 tissue_optimal_cut_off_max_gap_df, get_constraint = F, get_conserv = F,  
                                 propor_samples_split_read = 0.05, extend_region_to_plot = 10000, 
                                 collapseTranscripts, transcriptAnnotation, 
                                 aceview_annot = NULL, add_custom_annot_track = NULL, 
                                 all_split_reads = F){
  
  ##### convert symbol to ENSG ID #####
  
  if(str_detect(gene_id, "ENSG")){
    
    ensembl_gene_id <- gene_id
    
  }else{
    
    ensembl_gene_id_to_symbol_df_dup <- 
    ensembl_gene_id_to_symbol_df %>%
      filter(!duplicated(external_gene_name), # removing instances whereby 2 ensembl ids match onto 1 
             external_gene_name == gene_id)
  
    if(nrow(ensembl_gene_id_to_symbol_df_dup) == 0){
      
      stop(str_c(gene_id, " not found in the ensembl database (v92)"))
      
    }else{
      
      ensembl_gene_id <- ensembl_gene_id_to_symbol_df_dup$ensembl_gene_id
      
    }
    
  }
  
  ##### checking statements #####
  
  # extract all gene definitions from ensembl txdb
  ensembl_grch38_v92_genes_txdb_genes <- genes(ensembl_grch38_v92_genes_txdb)
  
  if(!ensembl_gene_id %in% ensembl_grch38_v92_genes_txdb_genes$gene_id){
    
    stop(str_c(ensembl_gene_id, " not found in the txdb provided"))
    
  }
  
  if(all(!tissues_to_plot %in% ERs_w_annotation_df$tissue)){
    
    stop(str_c(tissues_to_plot, " not found in tissues"))
    
  }
  
  ##### getting the start stop positions to plot #####
  
  gene_cord <- ensembl_grch38_v92_genes_txdb_genes[ensembl_grch38_v92_genes_txdb_genes$gene_id == ensembl_gene_id]
  seqnames_to_plot <- as.character(seqnames(gene_cord))
  start_to_plot <- (start(gene_cord) - extend_region_to_plot)
  end_to_plot <- (end(gene_cord) + extend_region_to_plot)
  bp_plotted <- end_to_plot - start_to_plot
  
  ERs_w_annotation_df_to_plot <- 
    ERs_w_annotation_df %>% 
    filter(seqnames == seqnames_to_plot, start >= start_to_plot, end <= end_to_plot, 
           tissue %in% tissues_to_plot)
  
  ##### Ideogram track - removed because not useful, only for aesthetics #####
  
  # ideo_track <- IdeogramTrack(genome = "hg38", chromosome = seqnames_to_plot, showBandId = TRUE, 
  #                             fontcolor = "black", size = 1)
  
  # plotTracks(ideo_track, from = start_to_plot, to = end_to_plot)
  
  ##### Genome axis track #####
  
  # get symbol to plot 
  gene_id_to_plot <- str_c(na.omit(c(ensembl_gene_id_to_symbol_df_dup$external_gene_name, ensembl_gene_id_to_symbol_df_dup$ensembl_gene_id)), collapse = "/")
  
  ga_track <- GenomeAxisTrack(range = IRanges(start = start(gene_cord), end = end(gene_cord), names = gene_id_to_plot), 
                              showId = TRUE, add53 = TRUE, add35 = TRUE, 
                              fill.range = get_palette("npg", 3)[3], col.range = "black", 
                              cex.id = 1, col = "black", fontcolor = "black", labelPos = "below", size = 2, 
                              background.title = "black", showTitle = F)
  
  # plotTracks(ga_track, from = start_to_plot, to = end_to_plot)
  
  ##### Annotation track - ERs #####
  
  annot_track_ERs_all_tissues <- get_annot_track_ERs(ERs_w_annotation_df_to_plot, tissues_to_plot)
  
  # plotTracks(annot_track_ERs_all_tissues[[1]], from = start_to_plot, to = end_to_plot)
  
  ##### Data track - split reads #####
  
  d_track_split_read_overlayed_all_tissues <- 
    get_data_track_split_read(tissues_to_plot = tissues_to_plot, gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df, 
                              ERs_w_annotation_df_to_plot = ERs_w_annotation_df_to_plot, gtex_split_read_table_annotated_only_junc_coverage = gtex_split_read_table_annotated_only_junc_coverage, 
                              propor_samples_split_read = propor_samples_split_read, gene_cord = gene_cord, extend_region_to_plot = extend_region_to_plot, all_split_reads = all_split_reads)
  
  ##### Merge ER and split read tracks #####
  
  ER_split_tracks_ordered_all_tissues <- merge_ER_split_read_tracks(tissues_to_plot, annot_track_ERs_all_tissues, d_track_split_read_overlayed_all_tissues)
  
  # plotTracks(ER_split_tracks_ordered_all_tissues, from = start_to_plot, to = end_to_plot, 
  #             sizes = rep(1, length(tissues_to_plot) * 4))
  
  ##### Gene region track - ensembl #####
  
  print("plotting gene region track for ensembl v92")
  
  gr_track <- GeneRegionTrack(ensembl_grch38_v92_genes_txdb, genome = genome_build, chromosome = seqnames_to_plot, name = "v92", 
                              collapseTranscripts = collapseTranscripts, transcriptAnnotation = transcriptAnnotation, from = start_to_plot, to = end_to_plot, 
                              background.title = "black", fontcolor.group = "black", col = "black", col.line = "black" , cex.group = 0.5, size = 1.75)
  
  all_annot_tracks <- c(ga_track, ER_split_tracks_ordered_all_tissues, gr_track)
  
  # plotTracks(gr_track, from = start_to_plot, to = end_to_plot, cex.group = 0.4)
  
  ##### Gene region track - aceview #####
  
  if(!is.null(aceview_annot)){
    
    if(class(aceview_annot) != "TxDb"){
      
      print("Warning: AceView input not a txdb, skipping adding this as a genome region track")
      
    }else{
      
      print("plotting gene region track for aceview")
      
      gr_track_aceview <- GeneRegionTrack(aceview_annot, genome = genome_build, chromosome = seqnames_to_plot, name = "AV", 
                                          collapseTranscripts = "meta", transcriptAnnotation = "gene", from = start_to_plot, to = end_to_plot, 
                                          background.title = "black", fontcolor.group = "black", col = "black", col.line = "black" , cex.group = 0.5, size = 1.75)
      
      all_annot_tracks <- c(all_annot_tracks, gr_track_aceview)
      
    }

    
  }
  
  # only plot mean coverage for single tissue
  if(length(tissues_to_plot) == 1){
    
    ##### get GR for data tracks - mean coverage/conservation #####
    
    range_to_plot_ir <- IRanges(start = start_to_plot:end_to_plot, width = 1)
    
    range_to_plot_gr <- GRanges(seqnames = seqnames_to_plot, 
                                range_to_plot_ir)
    
    ##### Data track - mean coverage #####
    
    d_track_mean_cov_w_cut_off <- get_data_track_mean_cov(gtex_split_read_table_mean_cov_df, tissues_to_plot, range_to_plot_ir, range_to_plot_gr, seqnames_to_plot)
    
    all_annot_tracks <- c(all_annot_tracks, d_track_mean_cov_w_cut_off)
    
  }
  
  
  ##### Data track - conservation #####
  
  if(get_conserv == T){
    
    range_to_plot_gr <-
      get_conservation_score_for_regions(conserv_score_to_load = "phast_cons_7", gr = range_to_plot_gr)
    
    d_track_conserv <- DataTrack(range_to_plot_gr, type = "a", aggregation = "mean", window = "fixed", windowSize = 25, background.title = "black",
                                 name = "phastcons", size = 1)
    
    # plotTracks(d_track_conserv, from = start_to_plot, to = end_to_plot)
    
    range_to_plot_gr$mean_phast_cons_7 <- NULL
    
    all_annot_tracks <- c(all_annot_tracks, d_track_conserv)
    
  }
  
  ##### Data track - constraint #####
  
  if(get_constraint == T){
    
    d_track_constraint <- get_data_track_constraint(seqnames_to_plot, start_to_plot, end_to_plot, CDTS_percentile_N7794_unrelated_all_chrs_gr)
    
    # plotTracks(d_track_constraint, from = start_to_plot, to = end_to_plot) 
    
    all_annot_tracks <- c(all_annot_tracks, d_track_constraint)
    
  }
  
  ##### Data track - custom #####
  
  if(!is.null(add_custom_annot_track)){
    
    annot_track_custom <- 
      AnnotationTrack(add_custom_annot_track, 
                      shape = "box", name = "custom annot", background.title = "black", size = 1, col = "black", fill = "black")
    
    # plotTracks(annot_track_custom, from = start_to_plot, to = end_to_plot, col = "black", fill = "black") 
    
    all_annot_tracks <- c(all_annot_tracks, annot_track_custom)
    
  }
  
  ##### plot tracks #####
  
  list_dummy_tracks <- get_dummy_track_titles(all_annot_tracks, gtex_split_read_table_mean_cov_df)
  
  plotTracks(trackList = all_annot_tracks, from = start_to_plot, to = end_to_plot, title.width = 0.4, fontsize = 8, 
             sizes = all_annot_tracks %>% lapply(FUN = function(x){ displayPars(x)$size}) %>% unlist())
  
  plotTracks(list_dummy_tracks, from = start_to_plot, to = end_to_plot, title.width = 0.3, add = T, fontsize = 12,
             sizes = list_dummy_tracks %>% lapply(FUN = function(x){ displayPars(x)$size}) %>% unlist())
            
}

##### Second level #####

source("/home/dzhang/projects/constraint_conservation_wd/constraint_conservation/conservation/conservation_general/conservation_general_functions.R")
source("/home/dzhang/misc_projects/bioinformatics_misc/bioinformatics_misc_git/query_biomart.R")

get_annot_track_ERs <- function(ERs_w_annotation_df_to_plot, tissues_to_plot){
  
  print(str_c("plotting ERs for ", str_c(tissues_to_plot, collapse = ", ")))
  
  annot_track_ERs_all_tissues <- list()
  
  for(i in seq_along(tissues_to_plot)){
    
    tissue_to_plot <- tissues_to_plot[i]
    
    ERs_w_annotation_gr <- 
      ERs_w_annotation_df_to_plot %>% 
      filter(tissue == tissue_to_plot) %>% 
      dplyr::select(seqnames, start, end, value, ensembl_grch38_v92_region_annot) %>% 
      as.data.frame() %>% 
      toGRanges()
    
    ERs_w_annotation_gr_annotated <- ERs_w_annotation_gr[str_detect(ERs_w_annotation_gr$ensembl_grch38_v92_region_annot, "exon")]
    ERs_w_annotation_gr_unannotated <- ERs_w_annotation_gr[ERs_w_annotation_gr$ensembl_grch38_v92_region_annot %in% c("intron", "intergenic")]
    
    annot_track_annotated_ERs <- AnnotationTrack(ERs_w_annotation_gr_annotated, 
                                                 shape = "box", name = tissue_to_plot, 
                                                 exon = "#003C67FF", `exon, intron` = "#0073C2FF", exon, intergenic = "#7AA6DCFF", 
                                                 intergenic = "#A73030FF", intron = "#CD534CFF", background.title = "black", 
                                                 col.line = "black", size = 1)
    
    annot_track_unannotated_ERs <- AnnotationTrack(ERs_w_annotation_gr_unannotated, 
                                                   shape = "box", name = tissue_to_plot,
                                                   exon = "#003C67FF", `exon, intron` = "#0073C2FF", exon, intergenic = "#7AA6DCFF", 
                                                   intergenic = "#A73030FF", intron = "#CD534CFF", stacking = "dense", background.title = "black", 
                                                   col.line = "black", size = 1)
    
    feature(annot_track_annotated_ERs) <-  ERs_w_annotation_gr_annotated$ensembl_grch38_v92_region_annot %>% as.character()                                             
    feature(annot_track_unannotated_ERs) <-  ERs_w_annotation_gr_unannotated$ensembl_grch38_v92_region_annot %>% as.character()
    
    annot_track_ERs_per_tissue <- list(annot_track_annotated_ERs, annot_track_unannotated_ERs)
    
    names(annot_track_ERs_per_tissue) <- c("annot", "unannot")
    
    annot_track_ERs_all_tissues[[i]] <- annot_track_ERs_per_tissue
    
    # plotTracks(trackList = annot_track_ERs_all_tissues[[1]], from = start_to_plot, to = end_to_plot)
    
  }
  
  names(annot_track_ERs_all_tissues) <- tissues_to_plot
  
  return(annot_track_ERs_all_tissues)
  
}

get_data_track_split_read <- function(tissues_to_plot, gtex_split_read_table_mean_cov_df, ERs_w_annotation_df_to_plot, 
                                      gtex_split_read_table_annotated_only_junc_coverage, 
                                      propor_samples_split_read, gene_cord, extend_region_to_plot, all_split_reads){
  
  print(str_c("plotting split reads for ", str_c(tissues_to_plot, collapse = ", ")))
  
  d_track_split_read_overlayed_all_tissues <- list()
  
  for(i in seq_along(tissues_to_plot)){
    
    tissue_to_plot <- tissues_to_plot[i]
    
    print(str_c(Sys.time(), " - plotting split reads for ", tissue_to_plot))
    
    gtex_split_read_table_annotated_path <- 
      gtex_split_read_table_mean_cov_df %>% 
      filter(OMIM_gtex_name == tissue_to_plot) %>% 
      .[["gtex_split_read_table_annotated_paths"]]
    
    tissue_n <- 
      gtex_split_read_table_mean_cov_df %>% 
      filter(OMIM_gtex_name == tissue_to_plot) %>% 
      .[["n_all"]]
    
    load(gtex_split_read_table_annotated_path)
    
    ERs_w_annotation_df_to_plot_tissue_filtered <- ERs_w_annotation_df_to_plot %>% filter(tissue == tissue_to_plot)
    
    d_track_split_read_overlayed <- 
      get_data_track_split_read_per_tissue(ERs_w_annotation_df_to_plot_tissue_filtered, tissue_to_plot, 
                                           gtex_split_read_table_annotated_only_junc_coverage, 
                                           propor_samples_split_read, tissue_n, gene_cord, extend_region_to_plot, all_split_reads)
    
    d_track_split_read_overlayed_all_tissues[[i]] <- d_track_split_read_overlayed
    
  }
  
  names(d_track_split_read_overlayed_all_tissues) <- tissues_to_plot
  
  return(d_track_split_read_overlayed_all_tissues)
  
}

merge_ER_split_read_tracks <- function(tissues_to_plot, annot_track_ERs_all_tissues, d_track_split_read_overlayed_all_tissues){
  
  stopifnot(identical(names(annot_track_ERs_all_tissues), names(d_track_split_read_overlayed_all_tissues)))
  
  ER_split_tracks_ordered_all_tissues <- list()
  
  for(i in seq_along(tissues_to_plot)){
    
    tissue_to_plot <- tissues_to_plot[i]
    
    ER_split_tracks_ordered_per_tissue <- 
      list(d_track_split_read_overlayed_all_tissues[[tissue_to_plot]], 
           annot_track_ERs_all_tissues[[tissue_to_plot]][["annot"]], annot_track_ERs_all_tissues[[tissue_to_plot]][["unannot"]])
           
    
    ER_split_tracks_ordered_all_tissues <- c(ER_split_tracks_ordered_all_tissues, ER_split_tracks_ordered_per_tissue)
    
  }
  
  ER_split_tracks_ordered_all_tissues_no_null <- ER_split_tracks_ordered_all_tissues[!(lapply(ER_split_tracks_ordered_all_tissues, is.null) %>% unlist())]
  
  return(ER_split_tracks_ordered_all_tissues_no_null)
  
}

get_data_track_mean_cov <- function(gtex_split_read_table_mean_cov_df, tissues_to_plot, range_to_plot_ir, range_to_plot_gr, seqnames_to_plot){
  
  gtex_mean_cov_path <- 
    gtex_split_read_table_mean_cov_df %>% 
    filter(OMIM_gtex_name == tissues_to_plot) %>% 
    .[["gtex_mean_cov_paths"]]
  
  gtex_mean_cov_tissue <- 
    gtex_mean_cov_path %>% str_replace("/.*/", "")
  
  gtex_mean_cov_path_chr <- 
    gtex_mean_cov_path %>% 
    str_c(., "/gtex_", gtex_mean_cov_tissue, "_", seqnames_to_plot, "_mean_cov.rda")
  
  load(gtex_mean_cov_path_chr)
  
  range_to_plot_gr$meanCoverage <- log(tissue_coverage_w_mean_normalised$meanCoverage[range_to_plot_ir], base = 10)
  
  d_track_mean_cov <- DataTrack(range_to_plot_gr, type = "a", aggregation = "mean", window = "fixed", windowSize = 25)
  
  # plotTracks(d_track_mean_cov, from = start_to_plot, to = end_to_plot)
  
  range_to_plot_gr$meanCoverage <- NULL
  
  tissue_optimal_cut_off_max_gap_df_cutoff <- 
    tissue_optimal_cut_off_max_gap_df %>% 
    filter(tissue == tissues_to_plot) %>% 
    .[["cut_off"]] %>% 
    log(base = 10)
  
  range_to_plot_gr$cut_off <- tissue_optimal_cut_off_max_gap_df_cutoff
  
  d_track_mean_cov_cut_off <- DataTrack(range_to_plot_gr[c(1, length(range_to_plot_gr))], type = "a", lty = 2, col = "red")
  
  d_track_mean_cov_w_cut_off <- OverlayTrack(trackList = list(d_track_mean_cov, d_track_mean_cov_cut_off), name = "log10(Mean coverage)", 
                                             size = 1, background.title = "black")
  
  return(d_track_mean_cov_w_cut_off)
  
}

get_data_track_constraint <- function(seqnames_to_plot, start_to_plot, end_to_plot, CDTS_percentile_N7794_unrelated_all_chrs_gr){
  
  CDTS_percentile_N7794_unrelated_all_chrs_gr_chr_filtered <- 
    CDTS_percentile_N7794_unrelated_all_chrs_gr %>% 
    keepSeqlevels(seqnames_to_plot, pruning.mode = "coarse")
  
  constraint_hits <- 
    findOverlaps(CDTS_percentile_N7794_unrelated_all_chrs_gr_chr_filtered, 
                 GRanges(seqnames = seqnames_to_plot, IRanges(start = start_to_plot, end = end_to_plot)), 
                 type = "within")
  
  CDTS_percentile_N7794_unrelated_all_chrs_gr_chr_region_filtered <- 
    CDTS_percentile_N7794_unrelated_all_chrs_gr_chr_filtered[queryHits(constraint_hits)]
  
  CDTS_percentile_N7794_unrelated_all_chrs_gr_chr_region_filtered$percentile <- NULL
  
  d_track_constraint <- DataTrack(CDTS_percentile_N7794_unrelated_all_chrs_gr_chr_region_filtered, 
                                  type = "a", aggregation = "mean", window = "fixed", windowSize = 25, background.title = "black",
                                  name = "CDTS", size = 1)
  
  # plotTracks(d_track_constraint, from = start_to_plot, to = end_to_plot) 
  
  return(d_track_constraint)
  
}

get_dummy_track_titles <- function(all_annot_tracks, gtex_split_read_table_mean_cov_df){
  
  name_size_track_class_list <- 
  all_annot_tracks %>% lapply(FUN = function(x){
    
    y <- x %>% displayPars()
    
    name_size_track_df <- data_frame(track_name = names(x), 
                                     track_size = y[[c("size")]], 
                                     track_type = class(x))
    
    return(name_size_track_df)
    
  })
  
  name_size_track_df <- 
    do.call("bind_rows", name_size_track_class_list) %>% 
    left_join(gtex_split_read_table_mean_cov_df %>% dplyr::select(OMIM_gtex_name, gtex_tissues_name_to_plot), by = c("track_name" = "OMIM_gtex_name")) %>% 
    mutate(track_name = ifelse(is.na(gtex_tissues_name_to_plot), track_name, gtex_tissues_name_to_plot)) %>% 
    group_by(track_name) %>% 
    mutate(total_track_size = sum(track_size)) %>% 
    filter(!duplicated(track_name)) 
  
  list_dummy_tracks <- list()
  
  for(i in 1:nrow(name_size_track_df)){
  
    track_name_to_plot <- name_size_track_df$track_name[i] %>% as.character()
    track_size_to_plot <- name_size_track_df$total_track_size[i]
    track_type_to_plot <- name_size_track_df$track_type[i]
      
    list_dummy_tracks[[i]] <- 
      DataTrack(GRanges("chrY:69-69"), genome = "hg38", 
                size = track_size_to_plot, 
                showAxis = F, 
                showTitle = T, 
                name = track_name_to_plot, 
                background.title = "black")

    displayPars(list_dummy_tracks[[i]])
    
  }
  
  return(list_dummy_tracks)
  
}

##### Third level #####

get_data_track_split_read_per_tissue <- function(ERs_w_annotation_df_to_plot_tissue_filtered, tissue_to_plot, 
                                                 gtex_split_read_table_annotated_only_junc_coverage, 
                                                 propor_samples_split_read, tissue_n, gene_cord, extend_region_to_plot, all_split_reads){
  
  annot_ERs <- 
    ERs_w_annotation_df_to_plot_tissue_filtered %>% 
    filter(str_detect(ensembl_grch38_v92_region_annot, "exon"))
  
  unannot_ERs <- 
    ERs_w_annotation_df_to_plot_tissue_filtered %>% 
    filter(ensembl_grch38_v92_region_annot %in% c("intron", "intergenic")) 
  
  annot_ERs_jun_ids <- 
    c((annot_ERs$p_annot_junc_ids_split_read_annot %>% str_split(";") %>% unlist()), 
      (annot_ERs$unannot_junc_ids_split_read_annot %>% str_split(";") %>% unlist())) %>% 
    unique() %>% 
    na.omit()
  
  unannot_ERs_jun_ids <- 
    c((unannot_ERs$p_annot_junc_ids_split_read_annot %>% str_split(";") %>% unlist()), 
      (unannot_ERs$unannot_junc_ids_split_read_annot %>% str_split(";") %>% unlist())) %>% 
    unique() %>% 
    na.omit()
  
  p_annot_junc_ids <- ERs_w_annotation_df_to_plot_tissue_filtered$p_annot_junc_ids_split_read_annot %>% str_split(";") %>% unlist() %>% unique()
  p_annot_junc_ids <- p_annot_junc_ids[!is.na(p_annot_junc_ids)]
  unannot_junc_ids <- ERs_w_annotation_df_to_plot_tissue_filtered$unannot_junc_ids_split_read_annot %>% str_split(";") %>% unlist() %>% unique()
  unannot_junc_ids <- unannot_junc_ids[!is.na(unannot_junc_ids)]
  
  all_junc_ids <- c(p_annot_junc_ids, unannot_junc_ids)
  
  stopifnot(all(c(annot_ERs_jun_ids, unannot_ERs_jun_ids) %in% all_junc_ids))
  
  if(all_split_reads == T){
    
    gtex_split_read_table_annotated_only_junc_coverage_to_plot <- 
      gtex_split_read_table_annotated_only_junc_coverage %>% 
      as_tibble() %>% 
      filter(chr == seqnames(gene_cord) %>% str_replace("chr", ""),
             start >= (start(gene_cord) - extend_region_to_plot), 
             stop <= (end(gene_cord) + extend_region_to_plot), 
             strand == strand(gene_cord) %>% as.character())
    
  }else{
    
    gtex_split_read_table_annotated_only_junc_coverage_to_plot <- 
      gtex_split_read_table_annotated_only_junc_coverage %>% 
      filter(junID %in% all_junc_ids, 
             strand == strand(gene_cord) %>% as.character())
    
  }
  
  if(nrow(gtex_split_read_table_annotated_only_junc_coverage_to_plot) == 0){
    
    print(str_c(Sys.time(), " - no split reads to plot attached to ERs in the ", tissue_to_plot))
    
    return(NULL)
    
  }
  
  gtex_split_read_table_annotated_only_junc_coverage_to_plot <-
    gtex_split_read_table_annotated_only_junc_coverage_to_plot %>%
    mutate(annot_unannot_ER = ifelse(junID %in% annot_ERs_jun_ids,
                                     ifelse(junID %in% unannot_ERs_jun_ids, "both", "annot"),
                                     ifelse(junID %in% unannot_ERs_jun_ids, "unannot", "none")),
           annot_unannot_split_read = ifelse(precBoundDonor == T | precBoundAcceptor == T, "p_annot", "unannot"),
           propor_samples = countsSamples/tissue_n) %>%
    filter(!(annot_unannot_ER == "none" & precBoundDonor == T & precBoundAcceptor == T))
  
  d_track_split_read_all <- list()
  
  for(j in 1:nrow(gtex_split_read_table_annotated_only_junc_coverage_to_plot)){
    
    split_read_start <- gtex_split_read_table_annotated_only_junc_coverage_to_plot$start[j]
    split_read_end <- gtex_split_read_table_annotated_only_junc_coverage_to_plot$stop[j]
    split_read_annot_uannot_ER <- gtex_split_read_table_annotated_only_junc_coverage_to_plot$annot_unannot_ER[j]
    split_read_annot_uannot_split_read <- gtex_split_read_table_annotated_only_junc_coverage_to_plot$annot_unannot_split_read[j]
    split_read_propor_samples <- gtex_split_read_table_annotated_only_junc_coverage_to_plot$propor_samples[j]
    
    split_read_gr_to_plot <- 
      GRanges(seqnames = str_c("chr", gtex_split_read_table_annotated_only_junc_coverage_to_plot$chr[j]), 
              IRanges(start = c(split_read_start, round((split_read_end+split_read_start)/2, digits = 0), split_read_end), width = 1))
    
    split_read_gr_to_plot$value <- c(0, 1, 0)
    
    if(split_read_annot_uannot_ER == "unannot"){
      
      line_colour <- "#A73030FF"
      
    }else if(split_read_annot_uannot_ER == "annot"){
      
      line_colour <- "#003C67FF"
      
    }else if(split_read_annot_uannot_ER == "both"){
      
      line_colour <- "#0F8C87"
      
    }else if(split_read_annot_uannot_ER == "none"){
      
      line_colour <- "black"
      
    }
    
    line_type <- ifelse(split_read_annot_uannot_split_read == "p_annot", 1, 2)
    
    d_track_split_read <- DataTrack(split_read_gr_to_plot, type = "a", lty = line_type, col = line_colour, showAxis = F, lwd = (5 * split_read_propor_samples))
    
    d_track_split_read_all[[j]] <- d_track_split_read
    
  }
  
  d_track_split_read_overlayed <- OverlayTrack(trackList = d_track_split_read_all, background.title = "black", size = 1.75, name = tissue_to_plot)
  
  
  return(d_track_split_read_overlayed)
  
}

# Main ------------------------------------------------------------------------------------------------

gtex_split_read_table_mean_cov_df <- get_gtex_split_read_table_mean_cov_n_samples_df(gtex_tissue_name_formatting)

ensembl_grch38_v92_genes_txdb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/ensembl/gtf_gff3/v92/Homo_sapiens.GRCh38.92.gtf", 
                         output_path = "/data/references/ensembl/txdb_sqlite/v92/ensembl_grch38_v92_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), genome_build = "hg38")

aceview_hg38_txdb <- 
  generate_txDb_from_gtf(gtf_gff3_path = "/data/references/aceview/gff/AceView.ncbi_37.genes_gff.gff", 
                         output_path = "/data/references/aceview/txdb_sqlite/aceview_hg38_txdb.sqlite",
                         seq_levels_to_keep = c(1:22, "X", "Y"), genome_build = "hg19", convert_hg19_to_hg38 = T)

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific <-
  ERs_w_annotation_all_tissues %>%
  filter(width > 3, !str_detect(tissue, "cells|testis|vagina|ovary|uterus|prostate|cervix|bladder|fallopian|breast"),
         ensembl_grch38_v92_region_annot != "exon, intergenic, intron")

# ERs_w_annotation_df <- ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific
# txdb <- ensembl_grch38_v92_genes_txdb
# ensembl_gene_id_to_symbol_df <- ensembl_gene_id_to_symbol_df_v92
# gene_id <- "SNCA"
# tissues_to_plot <- c("frontalcortexba9", "brain_cerebellum")
# genome_build <- "hg38"
# get_constraint <- T
# get_conserv <- F
# propor_samples_split_read <- 0.05
# extend_region_to_plot <- 30000
# aceview_annot <- aceview_hg38_txdb
# add_custom_annot_track <- NULL
# all_split_reads <- F
# 
# visualise_ER_example(ERs_w_annotation_df, txdb, ensembl_gene_id_to_symbol_df,
#                      ensembl_gene_id,
#                      tissues_to_plot, genome_build = "hg38",
#                      gtex_split_read_table_mean_cov_df,
#                      tissue_optimal_cut_off_max_gap_df,
#                      get_constraint,
#                      get_conserv,
#                      propor_samples_split_read,
#                      extend_region_to_plot,
#                      collapseTranscripts = "meta",
#                      transcriptAnnotation = "gene",
#                      aceview_annot = NULL,
#                      add_custom_annot_track = NULL,
#                      all_split_reads)
