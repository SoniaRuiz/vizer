library(tidyverse)
library(stringr)
library(Gviz)
library(ggpubr)
library(regioneR)
library(glue)

# Set WD ----------------------------------------------------------------------------------------------

if(!web_app){ 
  
  OMIM_wd <- Sys.getenv("OMIM_wd") 
  setwd(OMIM_wd)
  
}

# Load data -------------------------------------------------------------------------------------------

# need to be changed by sonia 
mean_coverage_path <- "/root/vizER/vizER_data/by_tissue_smfrze_use_me/"
path_to_data_folder <- "/root/vizER/vizER_data/"
output_image_path <- "/root/vizER/www/OMIM_reannot_plot.png"

constraint_grs_split_by_chr_paths_df <- 
  data_frame(paths = list.files(str_c(path_to_data_folder, "constraint_grs_split_by_chr/"), full.names = T), 
             chr = paths %>% str_replace(".*CDTS_percentile_N7794_unrelated_", "") %>% 
               str_replace("\\.rda", ""))

ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db <- 
  dbConnect(RSQLite::SQLite(), 
            str_c(path_to_data_folder, "ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific.sqlite"))

tissue_optimal_cut_off_max_gap_df <- read_delim(str_c(path_to_data_folder, "exon_delta_details_optimised_maxgap_cutoff.csv"), delim = ",")

gtex_tissue_name_formatting <- read_delim(str_c(path_to_data_folder, "OMIM_gtex_tissue_name_formatting.csv"), delim = ",")

ensembl_gene_id_to_symbol_df_v92 <- read_delim(str_c(path_to_data_folder, "ensembl_gene_id_to_symbol_df_v92.csv"), delim = ",")

# Functions -------------------------------------------------------------------------------------------

##### First level #####

source(str_c(path_to_data_folder, "generate_txDb_from_gtf.R"))

get_gtex_split_read_table_mean_cov_n_samples_df <- function(gtex_tissue_name_formatting){
  
  gtex_split_read_table_annotated_paths <- 
    list.files(str_c(path_to_data_folder, "gtex_split_read_table_annotated_rda/"), full.names = T)
  
  gtex_split_read_table_df <- 
    data_frame(gtex_split_read_table_annotated_paths = gtex_split_read_table_annotated_paths, 
               tissue = 
                 gtex_split_read_table_annotated_paths %>% 
                 str_replace(".*/", "") %>% 
                 str_replace("_split_read_table_annotated.rda", "") %>% 
                 str_replace("brain", "") %>% 
                 str_replace_all("_", ""))
  
  gtex_mean_cov_df <- 
    data_frame(gtex_mean_cov_paths = 
                 list.files(mean_coverage_path, full.names = T), 
               tissue = 
                 list.files(mean_coverage_path, full.names = F) %>% 
                 str_replace("brain", "") %>% 
                 str_replace_all("_", ""))

  gtex_split_read_table_mean_cov_df <- 
    gtex_split_read_table_df %>% 
    inner_join(gtex_mean_cov_df) %>% 
    left_join(gtex_tissue_name_formatting, by = c("tissue" = "gtex_tissue_name_simplified"))
  
  return(gtex_split_read_table_mean_cov_df)
  
}

#' @param ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db sqlite db with ER coordinates and annotation
#' @param txdb txdb object with the definitions of genes inside 
#' @param ensembl_gene_id_to_symbol_df df with all the ensembl gene ids matched to symbols using biomart
#' @param gene_id gene id to plot 
#' @param tissues_to_plot tissues to plot
#' @param gtex_split_read_table_mean_cov_df df that maps the tissue to it's split read or mean coverage path
get_ER_table_to_display <- function(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, 
                                    txdb, ensembl_gene_id_to_symbol_df, gene_id, 
                                    tissues_to_plot, gtex_split_read_table_mean_cov_df, extend_region_to_plot){
  
  ##### convert symbol to ENSG ID #####
  
  if(str_detect(gene_id, "ENSG")){
    
    ensembl_gene_id <- gene_id
    
    ensembl_gene_id_to_symbol_df_dup <- 
      ensembl_gene_id_to_symbol_df %>%
      filter(!duplicated(external_gene_name), # removing instances whereby 2 ensembl ids match onto 1 
             ensembl_gene_id == gene_id)
    
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
  
  if(extend_region_to_plot == "auto"){
    
    extend_region_to_plot_auto <- width(gene_cord)/10
    start_to_plot <- (start(gene_cord) - extend_region_to_plot_auto)
    end_to_plot <- (end(gene_cord) + extend_region_to_plot_auto)
    bp_plotted <- end_to_plot - start_to_plot
    
  }else{
    
    start_to_plot <- (start(gene_cord) - extend_region_to_plot)
    end_to_plot <- (end(gene_cord) + extend_region_to_plot)
    bp_plotted <- end_to_plot - start_to_plot
    
  }
  
  filtered_ERs_sql <- glue_sql('SELECT * FROM ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific WHERE "seqnames" == ({seqname*})  AND "start" >= ({start*}) AND "end" <= ({end*}) AND "tissue" IN ({tissues*})',
                               seqname = seqnames_to_plot,
                               start = start_to_plot,
                               end = end_to_plot,
                               tissues = tissues_to_plot,
                               .con = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db)
  
  suppressWarnings(
    filtered_ERs_queried <- dbSendQuery(conn = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, statement = filtered_ERs_sql)
  )
  
  ERs_w_annotation_df_to_display <- dbFetch(filtered_ERs_queried) %>% as_tibble()
  
  ERs_w_annotation_df_to_display <- 
    ERs_w_annotation_df_to_display %>% 
    mutate(misannot_type = ifelse(!is.na(uniq_genes_split_read_annot), "split_read",
                                  ifelse(!is.na(overlap_any_gene_v92_name), "overlap", 
                                         ifelse(nearest_any_gene_v92_distance <= 10000, "within_10Kb", "none"))), 
           associated_gene = ifelse(misannot_type == "split_read", uniq_genes_split_read_annot, 
                                    ifelse(misannot_type == "overlap", overlap_any_gene_v92_name, 
                                           ifelse(misannot_type == "within_10Kb", nearest_any_gene_v92_name, NA))))
  
  ERs_w_annotation_df_to_display_w_split_read_data <- data_frame()
  
  for(i in seq_along(tissues_to_plot)){
    
    tissue_to_plot <- tissues_to_plot[i]
    
    ERs_w_annotation_df_to_display_tissue_filtered <- ERs_w_annotation_df_to_display %>% filter(tissue == tissue_to_plot)
    
    load(gtex_split_read_table_mean_cov_df %>% 
         filter(OMIM_gtex_name == tissue_to_plot) %>% 
         .[["gtex_split_read_table_annotated_paths"]])
    
    tissue_n <- 
      gtex_split_read_table_mean_cov_df %>% 
      filter(OMIM_gtex_name == tissue_to_plot) %>% 
      .[["n_all"]]
    
    ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details <- 
      get_split_read_info_for_ER_details_df(ERs_w_annotation_df_to_display_tissue_filtered, tissue_n, gtex_split_read_table_annotated_only_junc_coverage)
      
    ERs_w_annotation_df_to_display_w_split_read_data <- 
      ERs_w_annotation_df_to_display_w_split_read_data %>% 
      bind_rows(ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details)
    
  }
  
  ERs_w_annotation_df_to_display_w_split_read_data <- 
  ERs_w_annotation_df_to_display_w_split_read_data %>% 
  dplyr::select(ER_chr = seqnames, ER_start = start, ER_end = end, ER_width = width, tissue, mean_coverage = value, ensembl_grch38_v92_region_annot, misannot_type, associated_gene,
                overlap_any_gene_v92_name, nearest_any_gene_v92_name, nearest_any_gene_v92_distance, 
                split_read_annotation_type = annotationType_split_read_annot, annotationType_split_read_annot, split_read_to_any_gene = uniq_genes_split_read_annot,
                split_read_count_samp, split_read_propor_samp, split_read_starts, split_read_ends,
                mean_CDTS_percentile, mean_phast_cons_7, mean_phast_cons_100)
  
  return(ERs_w_annotation_df_to_display_w_split_read_data)
  
}

#' @param ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db sqlite db with ER coordinates and annotation 
#' @param txdb txdb object with the definitions of genes inside 
#' @param ensembl_gene_id_to_symbol_df df with all the ensembl gene ids matched to symbols using biomart
#' @param gene_id gene id to plot 
#' @param tissues_to_plot tissues to plot
#' @param genome_build genome build - defaults to hg38 
#' @param gtex_split_read_table_mean_cov_df df that maps the tissue to it's split read or mean coverage path
#' @param get_constraint lgl vector indicating whether to plot constraint scores
#' @param get_conserv lgl vector indicating whether to plot conserv scores 
#' @param get_mean_cov lgl vector indicating whether to plot mean cov (only works for 1 tissue)
#' @param propor_samples_split_read percentage of samples required that have the split read to plot
#' @param extend_region_to_plot number of bps either side of the gene to add to around plot
#' @param aceview_annot defaults to NULL otherwise, takes the aceview txdb 
#' @param add_custom_annot_track here gives the option to input your own positions to plot as a GRanges
#' @param all_split_reads lgl vector determining whether to plot all split reads in the region, not just those that connect to ERs
#' @return Gviz plot to visualise the ERs in the genomic region or by gene
visualise_ER_example <- function(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, txdb, ensembl_gene_id_to_symbol_df, gene_id, tissues_to_plot, 
                                 genome_build = "hg38", gtex_split_read_table_mean_cov_df, 
                                 tissue_optimal_cut_off_max_gap_df, get_constraint = F, get_conserv = F, get_mean_cov = F, 
                                 propor_samples_split_read = 0.05, extend_region_to_plot = "auto", 
                                 collapseTranscripts, transcriptAnnotation, 
                                 aceview_annot = NULL, add_custom_annot_track = NULL, 
                                 all_split_reads = F){
  
  ##### convert symbol to ENSG ID #####
  
  if(str_detect(gene_id, "ENSG")){
    
    ensembl_gene_id <- gene_id
    
    ensembl_gene_id_to_symbol_df_dup <- 
      ensembl_gene_id_to_symbol_df %>%
      filter(!duplicated(external_gene_name), # removing instances whereby 2 ensembl ids match onto 1 
             ensembl_gene_id == gene_id)
    
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
  
  if(web_app){ setProgress(value = 0.14) }
  
  ##### checking statements #####
  
  # extract all gene definitions from ensembl txdb
  ensembl_grch38_v92_genes_txdb_genes <- genes(ensembl_grch38_v92_genes_txdb)
  
  if(!ensembl_gene_id %in% ensembl_grch38_v92_genes_txdb_genes$gene_id){
    
    stop(str_c(ensembl_gene_id, " not found in the txdb provided"))
    
  }
  
  ##### getting the start stop positions to plot #####
  
  gene_cord <- ensembl_grch38_v92_genes_txdb_genes[ensembl_grch38_v92_genes_txdb_genes$gene_id == ensembl_gene_id]
  seqnames_to_plot <- as.character(seqnames(gene_cord))
  
  if(extend_region_to_plot == "auto"){
    
    extend_region_to_plot_auto <- width(gene_cord)/10
    start_to_plot <- (start(gene_cord) - extend_region_to_plot_auto)
    end_to_plot <- (end(gene_cord) + extend_region_to_plot_auto)
    bp_plotted <- end_to_plot - start_to_plot
    
  }else{
    
    start_to_plot <- (start(gene_cord) - extend_region_to_plot)
    end_to_plot <- (end(gene_cord) + extend_region_to_plot)
    bp_plotted <- end_to_plot - start_to_plot
    
  }
  
  if(web_app){ setProgress(value = 0.18) }

  filtered_ERs_sql <- glue_sql('SELECT * FROM ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific WHERE "seqnames" == ({seqname*})  AND "start" >= ({start*}) AND "end" <= ({end*}) AND "tissue" IN ({tissues*})',
                               seqname = seqnames_to_plot,
                               start = start_to_plot,
                               end = end_to_plot,
                               tissues = tissues_to_plot,
                               .con = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db)
  
  suppressWarnings(
    filtered_ERs_queried <- dbSendQuery(conn = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, statement = filtered_ERs_sql)
  )
  
  ERs_w_annotation_df_to_plot <- dbFetch(filtered_ERs_queried) %>% as_tibble()
  
  if(web_app){ setProgress(value = 0.2) }
  
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
                              cex.id = 0.8, col = "black", fontcolor = "black", labelPos = "below", size = 2, 
                              background.title = "black", showTitle = T, name = seqnames_to_plot)
  
  if(web_app){ setProgress(value = 0.26) }
  
  # plotTracks(ga_track, from = start_to_plot, to = end_to_plot)
  
  ##### Annotation track - ERs #####
  
  annot_track_ERs_all_tissues <- get_annot_track_ERs(ERs_w_annotation_df_to_plot, tissues_to_plot)
  
  if(web_app){ setProgress(value = 0.35) }
  
  # plotTracks(annot_track_ERs_all_tissues[[1]], from = start_to_plot, to = end_to_plot)
  
  ##### Data track - split reads #####
  
  d_track_split_read_overlayed_all_tissues <- 
    get_data_track_split_read(tissues_to_plot = tissues_to_plot, gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df, 
                              ERs_w_annotation_df_to_plot = ERs_w_annotation_df_to_plot, gtex_split_read_table_annotated_only_junc_coverage = gtex_split_read_table_annotated_only_junc_coverage, 
                              propor_samples_split_read = propor_samples_split_read, gene_cord = gene_cord, extend_region_to_plot = extend_region_to_plot, all_split_reads = all_split_reads)
  
  if(web_app){ setProgress(value = 0.55) }
  
  ##### Merge ER and split read tracks #####
  
  ER_split_tracks_ordered_all_tissues <- merge_ER_split_read_tracks(tissues_to_plot, annot_track_ERs_all_tissues, d_track_split_read_overlayed_all_tissues)
  
  # plotTracks(ER_split_tracks_ordered_all_tissues, from = start_to_plot, to = end_to_plot, 
  #             sizes = rep(1, length(tissues_to_plot) * 4))
  
  ##### Gene region track - ensembl #####
  
  print("plotting gene region track for ensembl v92")
  
  gr_track <- GeneRegionTrack(ensembl_grch38_v92_genes_txdb, genome = genome_build, chromosome = seqnames_to_plot, name = "Ens92", 
                              collapseTranscripts = collapseTranscripts, transcriptAnnotation = transcriptAnnotation, from = start_to_plot, to = end_to_plot, 
                              background.title = "black", fontcolor.group = "black", col = "black", col.line = "black" , cex.group = 0.5, size = 1.75)
  
  if(web_app){ setProgress(value = 0.65) }
  
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
  
  range_to_plot_ir <- IRanges(start = start_to_plot:end_to_plot, width = 1)
  
  range_to_plot_gr <- GRanges(seqnames = seqnames_to_plot, 
                              range_to_plot_ir)
  
  # only plot mean coverage for single tissue
  if(length(tissues_to_plot) == 1){
    
    if(get_mean_cov == T){
      
      ##### Data track - mean coverage #####
      
      d_track_mean_cov_w_cut_off <- get_data_track_mean_cov(gtex_split_read_table_mean_cov_df, tissues_to_plot, range_to_plot_ir, range_to_plot_gr, seqnames_to_plot)
      
      all_annot_tracks <- c(all_annot_tracks, d_track_mean_cov_w_cut_off)
      
    }
    
  }
  
  
  ##### Data track - conservation #####
  
  if(get_conserv == T){
    
    suppressWarnings(
      range_to_plot_gr_conserv <-
        get_conservation_score_for_regions_bw(bw_path = str_c(path_to_data_folder, "hg38.phastCons7way.bw"), gr = range_to_plot_gr)
    )
    
    d_track_conserv <- DataTrack(range_to_plot_gr_conserv, type = "a", aggregation = "mean", window = "fixed", windowSize = 25, background.title = "black",
                                 name = "PC7", size = 1)
    
    # plotTracks(d_track_conserv, from = start_to_plot, to = end_to_plot)
    
    range_to_plot_gr_conserv$mean_phastCons7way <- NULL
    
    all_annot_tracks <- c(all_annot_tracks, d_track_conserv)
    
  }
  
  ##### Data track - constraint #####
  
  if(get_constraint == T){
    
    d_track_constraint <- get_data_track_constraint(seqnames_to_plot, start_to_plot, end_to_plot, constraint_grs_split_by_chr_paths_df)
    
    # plotTracks(d_track_constraint, from = start_to_plot, to = end_to_plot) 
    
    all_annot_tracks <- c(all_annot_tracks, d_track_constraint)
    
  }
  
  ##### Get dummy tracks #####
  
  list_dummy_tracks <- get_dummy_track_titles(all_annot_tracks, gtex_split_read_table_mean_cov_df)
  
  if(web_app){ setProgress(value = 0.85) }
  
  ##### Data track - custom #####
  
  if(!is.null(add_custom_annot_track) && add_custom_annot_track != ""){
    
    add_custom_annot_track_gr <- GRanges(add_custom_annot_track)
    
    add_custom_annot_track_overlaps <- findOverlaps(add_custom_annot_track_gr, GRanges(str_c(seqnames_to_plot, ":", 
                                                                                          start_to_plot, "-", 
                                                                                          end_to_plot)))
    
    if(length(queryHits(add_custom_annot_track_overlaps)) != length(add_custom_annot_track)) stop("some custom track SNPs do not overlap with region to plot")
    
    all_annot_tracks_highlighted <- HighlightTrack(trackList = all_annot_tracks, range = add_custom_annot_track_gr)

    if(web_app){ png(output_image_path, res = 600, width = 10, height = 11.69/2, units = "in") }
    
    plotTracks(trackList = all_annot_tracks_highlighted, from = start_to_plot, to = end_to_plot, title.width = 0.5, fontsize = 7, 
               sizes = all_annot_tracks %>% lapply(FUN = function(x){ displayPars(x)$size}) %>% unlist())
    
    plotTracks(list_dummy_tracks, from = start_to_plot, to = end_to_plot, title.width = 0.5, add = T, fontsize = 7,
               sizes = list_dummy_tracks %>% lapply(FUN = function(x){ displayPars(x)$size}) %>% unlist())
    if(web_app){ dev.off() }
    
  }else{
    
    if(web_app){ png(output_image_path, res = 600, width = 10, height = 11.69/2, units = "in") }
    plotTracks(trackList = all_annot_tracks, from = start_to_plot, to = end_to_plot, title.width = 0.5, fontsize = 7, 
               sizes = all_annot_tracks %>% lapply(FUN = function(x){ displayPars(x)$size}) %>% unlist())
    
    plotTracks(list_dummy_tracks, from = start_to_plot, to = end_to_plot, title.width = 0.5, add = T, fontsize = 7,
               sizes = list_dummy_tracks %>% lapply(FUN = function(x){ displayPars(x)$size}) %>% unlist())
    if(web_app){ dev.off() }
    
  }


            
}

##### Second level #####

source(str_c(path_to_data_folder, "conservation_general_functions_from_bw.R"))
source(str_c(path_to_data_folder, "query_biomart.R"))

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
  
  ER_split_tracks_ordered_all_tissues_no_na <- ER_split_tracks_ordered_all_tissues[!is.na(ER_split_tracks_ordered_all_tissues)]
  
  return(ER_split_tracks_ordered_all_tissues_no_na)
  
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
  
  d_track_mean_cov_cut_off <- DataTrack(range_to_plot_gr, type = "a", lty = 2, col = "red", type = "a", aggregation = "mean", window = "auto")
  
  d_track_mean_cov_w_cut_off <- OverlayTrack(trackList = list(d_track_mean_cov, d_track_mean_cov_cut_off), name = "MC", 
                                             size = 1, background.title = "black")
  
  return(d_track_mean_cov_w_cut_off)
  
}

get_data_track_constraint <- function(seqnames_to_plot, start_to_plot, end_to_plot, constraint_grs_split_by_chr_paths_df){
  
  load(constraint_grs_split_by_chr_paths_df %>% 
         filter(chr == seqnames_to_plot) %>% 
         .[["paths"]])
  
  constraint_hits <- 
    findOverlaps(CDTS_percentile_N7794_unrelated_one_chr_gr, 
                 GRanges(seqnames = seqnames_to_plot, IRanges(start = start_to_plot, end = end_to_plot)), 
                 type = "within")
  
  CDTS_percentile_N7794_unrelated_one_chr_gr_region_filtered <- 
    CDTS_percentile_N7794_unrelated_one_chr_gr[queryHits(constraint_hits)]
  
  CDTS_percentile_N7794_unrelated_one_chr_gr_region_filtered$percentile <- NULL
  
  d_track_constraint <- DataTrack(CDTS_percentile_N7794_unrelated_one_chr_gr_region_filtered, 
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

get_split_read_info_for_ER_details_df <- function(ERs_w_annotation_df_to_display_tissue_filtered, tissue_n, gtex_split_read_table_annotated_only_junc_coverage){
  
  ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details <- 
    ERs_w_annotation_df_to_display_tissue_filtered %>% mutate(split_read_count_samp = as.integer(NA), 
                                                              split_read_propor_samp = as.numeric(NA),
                                                              split_read_starts = as.integer(NA), 
                                                              split_read_ends = as.integer(NA))
  
  gtex_split_read_table_annotated_only_junc_coverage_int_jun_id <- 
      gtex_split_read_table_annotated_only_junc_coverage %>% as_tibble() %>% mutate(junID = junID %>% as.integer()) %>% arrange(junID)
  
  for(j in 1:nrow(ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details)){
    
    p_annot_junc_ids_split_reads_int <- ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$p_annot_junc_ids_split_read_annot[j] %>% 
      str_split(";") %>% unlist() %>% as.integer() %>% sort()
    
    if(length(p_annot_junc_ids_split_reads_int) == 0){
      
      next
      
    }else{
      
      gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered <- 
      gtex_split_read_table_annotated_only_junc_coverage_int_jun_id %>% 
        filter(junID %in% p_annot_junc_ids_split_reads_int)
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_count_samp[j] <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$countsSamples %>% as.integer() %>% str_c(collapse = ";")
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_propor_samp[j] <- 
        ((gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$countsSamples %>% as.integer())/tissue_n) %>% round(digits = 3) %>% str_c(collapse = ";")
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_starts[j] <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$start %>% str_c(collapse = ";")
      
      ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details$split_read_ends[j] <- 
        gtex_split_read_table_annotated_only_junc_coverage_int_jun_id_filtered$stop %>% str_c(collapse = ";")
      
    }
    
  }

  return(ERs_w_annotation_df_to_display_tissue_filtered_w_split_read_details)
    
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
    
    return(NA)
    
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
  generate_txDb_from_gtf(gtf_gff3_path = str_c(path_to_data_folder, "Homo_sapiens.GRCh38.92.gtf"), 
                         output_path = str_c(path_to_data_folder, "ensembl_grch38_v92_txdb.sqlite"),
                         seq_levels_to_keep = c(1:22, "X", "Y", "MT"), genome_build = "hg38")

aceview_hg38_txdb <- 
  generate_txDb_from_gtf(gtf_gff3_path = str_c(path_to_data_folder, "AceView.ncbi_37.genes_gff.gff"), 
                         output_path = str_c(path_to_data_folder, "aceview_hg38_txdb.sqlite"),
                         seq_levels_to_keep = c(1:22, "X", "Y"), genome_build = "hg19", convert_hg19_to_hg38 = T)

# ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db <- ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db
# txdb <- ensembl_grch38_v92_genes_txdb
# ensembl_gene_id_to_symbol_df <- ensembl_gene_id_to_symbol_df_v92
# gene_id <- "ERLIN1"
# tissues_to_plot <- c("brain_cerebellar_hemisphere")
# genome_build <- "hg38"
# get_constraint <- T
# get_conserv <- T
# get_mean_cov <- T
# propor_samples_split_read <- 0.05
# extend_region_to_plot <- "auto"
# collapseTranscripts <-  "meta"
# transcriptAnnotation <-  "gene"
# aceview_annot <- NULL
# add_custom_annot_track <- "chr10:100154922-100154922"
# all_split_reads <- F
# 
# ERs_w_annotation_df_to_display_w_split_read_data <-
#   get_ER_table_to_display(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, txdb, ensembl_gene_id_to_symbol_df, gene_id,
#                           tissues_to_plot, gtex_split_read_table_mean_cov_df, extend_region_to_plot)
# 
# visualise_ER_example(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db,
#                      txdb,
#                      ensembl_gene_id_to_symbol_df,
#                      gene_id,
#                      tissues_to_plot,
#                      genome_build,
#                      gtex_split_read_table_mean_cov_df,
#                      tissue_optimal_cut_off_max_gap_df,
#                      get_constraint,
#                      get_conserv,
#                      get_mean_cov,
#                      propor_samples_split_read,
#                      extend_region_to_plot,
#                      collapseTranscripts,
#                      transcriptAnnotation,
#                      aceview_annot,
#                      add_custom_annot_track,
#                      all_split_reads)
# 
# dev.print(file = "/home/dzhang/projects/OMIM_wd/OMIM_paper/web_application/ERLIN1_OMIM_reannot_example_vizER.png", device = png, res = 600, width = 10, height = 11.69/2, units = "in")

