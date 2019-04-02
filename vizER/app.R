#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyjs)
library(magrittr)
library(RSQLite)

# Generate alphabetical tissue choices

tissue_GTEx_choices <- c("Adipose - subcutaneous" =	"adipose_subcutaneous",
                         "Adipose - visceral" =	"adipose_visceral_omentum",
                         "Adrenal gland" =	"adrenal_gland",
                         "Aorta" =	"artery_aorta",
                         "Artery - coronary" =	"artery_coronary",
                         "Artery - tibial" =	"artery_tibial",
                         "Amygdala" =	"amygdala",
                         "Anterior cingulate cortex" =	"anteriorcingulatecortexba24",
                         "Caudate" =	"caudatebasalganglia",
                         "Cerebellar hemisphere" =	"brain_cerebellar_hemisphere",
                         "Cerebellum" =	"brain_cerebellum",
                         "Cortex"	= "cortex",
                         "Frontal Cortex" =	"frontalcortexba9",
                         "Hippocampus" =	"hippocampus",
                         "Hypothalamus"	= "hypothalamus",
                         "Nucleus accumbens" =	"nucleusaccumbensbasalganglia",
                         "Putamen" =	"putamenbasalganglia",
                         "Spinal cord" =	"spinalcordcervicalc-1",
                         "Substantia nigra" =	"substantianigra",
                         "Sigmoid" =	"colon_sigmoid",
                         "Transverse" =	"colon_transverse",
                         "Gastroesophageal junction" =	"esophagus_gastroesophageal_junction",
                         "Mucosa" =	"esophagus_mucosa",
                         "Muscularis" =	"esophagus_muscularis",
                         "Atrial appendage" =	"heart_atrial_appendage",
                         "Left ventricle" =	"heart_left_ventricle",
                         "Kidney" =	"kidney_cortex",
                         "Liver" =	"liver",
                         "Lung" =	"lung",
                         "Minor salivary gland" =	"minor_salivary_gland",
                         "Skeletal muscle" =	"muscle_skeletal",
                         "Nerve - tibial" =	"nerve_tibial",
                         "Pancreas" =	"pancreas",
                         "Pituitary" =	"pituitary",
                         "Skin (suprapubic)" =	"skin_not_sun_exposed_suprapubic",
                         "Skin (lower leg)"	= "skin_sun_exposed_lower_leg",
                         "Small Intestine" =	"small_intestine_terminal_ileum",
                         "Spleen" =	"spleen",
                         "Stomach" =	"stomach",
                         "Thyroid" =	"thyroid",
                         "Whole blood" =	"whole_blood")

tissue_GTEx_choices_alphabetical <- tissue_GTEx_choices[names(tissue_GTEx_choices) %>% order()]

# Set WD ----------------------------------------------------------------------------------------------

# OMIM_wd <- Sys.getenv("OMIM_wd")
 setwd(".")

# Functions -------------------------------------------------------------------------------------------

source("global.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  shiny::tags$head(
    includeScript("https://www.googletagmanager.com/gtag/js?id=UA-129116044-1"),
    shiny::tags$link(rel="shortcut icon", href="ucl-icon.png"),
    shiny::tags$script(src = "google-analytics.js"),

    shiny::tags$script('$( document ).ready(function() {
      /**********************************************************/
      /**** REDIRECT PARENT (in order to avoid iframe issue) ****/
      /**********************************************************/
      var path = $("#shinyframe", window.parent.document).attr("src");
    	if(path != undefined)
    		window.top.location.href = path;
    });'),
   



   shiny::tags$script('$(document).ready(function(){
  $(window).resize(function(){
    var $zoomImg = $("#plot img");
	var height = $zoomImg.height();

	$(".zoomWrapper").css("height", height);
	$(".zoomContainer .zoomWindow").css({"height": height});

	$.removeData($("#plot img"), "elevateZoom");
	$(".zoomContainer").remove();

	
	if (/Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini|Opera Mobile|Kindle|Windows Phone|PSP|AvantGo|Atomic Web Browser|Blazer|Puffin|QQbrowser|SEMC Browser|Skyfire|Tear|TeaShark|UC Browser|uZard Web|wOSBrowser|Yandex.Browser mobile/i.test(navigator.userAgent)) {
	$("#plot img").elevateZoom({
	zoomWindowPosition:7,
	zoomWindowWidth:height/2,
	zoomWindowHeight:height/2,
	responsive:true
	//zoomLensWidth:75,
	//zoomLensHeight:75
	});
	}
	else{
	$("#plot img").elevateZoom({
	zoomWindowPosition:11,
	zoomWindowWidth:height,
	zoomWindowHeight:height,
	scrollZoom:true,
	responsive:true
	});
	}
  });
});'),












    
    






































 
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  navbarPage(title = "Visualisation of Expressed Regions", 
             tabPanel(title = "vizER",
                      # Application title
                      titlePanel(title = "vizER"),
                      
                      sidebarLayout(
                        sidebarPanel(
                          textInput(inputId = "geneid", label = "Gene ID", value = "ERLIN1", width = NULL, placeholder = "Enter gene of interest (ENSG ID or gene symbol)"),
                          selectizeInput("tissue", "Tissue:", choices = tissue_GTEx_choices_alphabetical,  multiple = T, options = list(maxItems = 3), selected = "brain_cerebellar_hemisphere"),
                          h5(strong("Extend region to plot:")),
                          sliderInput(inputId = "extend_region_to_plot", label = checkboxInput("auto", "Auto", FALSE),min = 1000, max = 50000, value = 1000, step = 500),
                          checkboxInput(inputId = "get_mean_cov", label = "Plot mean coverage", value = FALSE),
                          checkboxInput(inputId = "get_conserv", label = "Plot conservation", value = FALSE),
                          checkboxInput(inputId = "get_constraint", label = "Plot constraint", value = FALSE),
                          textInput(inputId = "add_custom_annot_track", label = "SNP of interest:", value = "", 
                                    width = NULL, placeholder = "e.g. chr10:100154922-100154922"),
                          actionButton("update", "Accept")
                        ),
                        # Show a plot of the generated distribution
                        mainPanel(
                          tabsetPanel(
                            id = "gen_browser_panel",
                            selected = "Plot",
                            tabPanel(title = 'Plot',
                                     fillRow(
                                       column(width = 12, 
                                              shiny::tags$head(shiny::tags$script(src = "jquery.elevatezoom.js")),
                                              br(),
                                              singleton(
                                                shiny::tags$head(shiny::tags$script('Shiny.addCustomMessageHandler("startzoom",
                                                                                    function(message) {
                                                                                    if (/Android|webOS|iPhone|iPad|iPod|BlackBerry|IEMobile|Opera Mini|Opera Mobile|Kindle|Windows Phone|PSP|AvantGo|Atomic Web Browser|Blazer|Chrome Mobile|Dolphin|Dolfin|Doris|GO Browser|Jasmine|MicroB|Mobile Firefox|Mobile Safari|Mobile Silk|Motorola Internet Browser|NetFront|NineSky|Nokia Web Browser|Obigo|Openwave Mobile Browser|Palm Pre web browser|Polaris|PS Vita browser|Puffin|QQbrowser|SEMC Browser|Skyfire|Tear|TeaShark|UC Browser|uZard Web|wOSBrowser|Yandex.Browser mobile/i.test(navigator.userAgent)) { 
                                                                                    $("#plot img").elevateZoom({
                                                                                    zoomWindowPosition:7,
                                                                                    zoomWindowWidth:250,
                                                                                    zoomWindowHeight:200, 
										    responsive:true
                                                                                    //zoomLensWidth:75,
                                                                                    //zoomLensHeight:75
                                                                                    });
                                                                                    }
                                                                                    else{
                                                                                    $("#plot img").elevateZoom({
                                                                                    zoomWindowPosition:11,
                                                                                    zoomWindowWidth:450,
                                                                                    zoomWindowHeight:450,
                                                                                    scrollZoom:true,
										    responsive:true
                                                                                    });
                                                                                    }
                                                                                    }
                                                );'),
                                                 shiny::tags$script('Shiny.addCustomMessageHandler("stopzoom",
                                                                    function(message) {
                                                                    $.removeData($("#plot img"), "elevateZoom");
                                                                    $(".zoomContainer").remove();
                                                                    }
                                                 );'),
                                                 shiny::tags$script('$(document).on("shiny:value", function(e) {
                                                                    console.log(e.name);
                                                                    if (e.name == "plot") {  // mytable is the name / id of the output element
                                                                    console.log(e.name);
                                                                    Shiny.onInputChange("startzoom", "");
                                                                    }});')
                                                 )),
                                              
                                              imageOutput("plot")))),
                            tabPanel(title = 'ER summary table',
                                     fluidRow(
                                       column(width = 12, 
                                              br(),  
                                              DT::dataTableOutput("summary")))),
                            tabPanel(title = 'Download',
                                     fluidRow(
                                       column(width = 12, 
                                              br(),
                                              h3("Download page"),br(),
                                              h4("Plot"),
                                              p("Download plot in PNG format:"),
                                              downloadButton(outputId = "download_plot", label = "Download"), br(),br(),
                                              h4("ER summary table"),
                                              p("Download ER summary table in CSV format:"),
                                              downloadButton(outputId = "download_data", label = "Download"), br(),br(),
                                              h4("ERs from each tissue"),
                                              p("If you would like to obtain the ERs from each tissue please contact Mina Ryten (",a(href="mailto:mina.ryten@ucl.ac.uk","mina.ryten@ucl.ac.uk"),")."),
                                              p("We are working on making this data suitable for distribution. It will be ready soon."),
                                              downloadButton(outputId = "download_bed", label = "Download"))))

              )))),
              tabPanel("About",
                      titlePanel(title = "About"),
                      fluidPage(
                        navlistPanel(
                          tabPanel("Overview", fluidPage(
                            h3("Overview"),
                            p("Next-generation sequencing has equipped researchers to discover mutations across the entire genome. However, the knowledge required to decipher which of these mutations is causal for each patient’s disease is still incomplete at both gene and variant level, in particular for those mutations that fall into non-coding regions of the genome. We improve upon the existing annotation of majority of genes that are currently known to cause Mendelian disorders (OMIM genes) broadening the genetic horizon which can be used to prioritise variants and eventually, assign pathogenicity. With this in mind, we have developed this online web resource vizER with the primary goal of aiding clinical scientists and clinicians to visualise misannotations of any gene of interest, enabling better variant prioritisation and as a result, diagnosis of Mendelian disorders.")
                          )),
                          tabPanel("Publication", fluidPage(
                            h3("Abstract"),
                            p("Although the increasing use of whole-exome and whole-genome sequencing have improved the yield of genetic testing for Mendelian disorders, an estimated 50% of patients still leave the clinic without a genetic diagnosis. This can be is attributed in part to our lack of ability to accurately interpret the genetic variation detected through next-generation sequencing. Variant interpretation is fundamentally reliant on accurate and complete gene annotation, however numerous reports and discrepancies between gene annotation databases reveals that the knowledge of gene annotation remains far from comprehensive. Here, we detect and validate transcription in an annotation-agnostic manner across all 41 different GTEx tissues, then connect novel transcription to known genes, ultimately improving the annotation of 63% of the known OMIM-morbid genes. We find the majority of novel transcription to be tissue-specific in origin with brain tissues being most susceptible to misannotation. Furthermore, we find that novel transcribed regions tend to be poorly conserved, but are significantly depleted for genetic variation within humans suggesting they are functionally significant and potentially have human-specific functions. We present our findings through an online platform vizER, which enables individual genes to be visualised and queried for evidence of misannotation. We also release all tissue-specific transcriptomes in a BED format for ease of integration with whole-genome sequencing data. We anticipate that these resources will improve the diagnostic yield for a wide range of Mendelian disorders."), 
                            a(href="https://www.biorxiv.org/content/10.1101/499103v1", target="_blank","Read more")
                          ))#,
                          #tabPanel("Citation", fluidPage(
                          #  h3("Citation"),
                          #  p("To reference this resource please use:"), 
                          #  p("ref_for_OMIM_paper")
                          #))
              ))),
             tabPanel("Help",
                      fluidRow(
                        column(12,
                               h1("Help"),
                               p("Details of the parameters used as inputs and output plots and tables are described below:"),br(),
                               # h1("Annotated Properties"),
                               # tableOutput('annotated_properties'),br(),
                               h3("Input parameters"),
                               p("This is designed to allow you query individual genes of interest for misannotations from any of the 41 GTEx tissues and highlight a variant of interest if desired."),
                               shiny::tags$table(id = "help-input-data",
                                                 shiny::tags$tr(
                                                   shiny::tags$th("Input name"),
                                                   shiny::tags$th("Description")
                                                   ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Gene ID"),
                                                   shiny::tags$td("Gene symbol or Ensembl ID (based on version 92)")
                                                   ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Tissue"),
                                                   shiny::tags$td("Selection of which GTEx tissues to plot (maximum of 3 for each plot)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Extend region to plot"),
                                                   shiny::tags$td("Expands the region to plot defaulting to the +/- 10% of the length of the gene of interest (Auto) or a numeric value between 1000-50,000 indicating the number of bps")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Plot mean coverage"),
                                                   shiny::tags$td("Tick to plot the base-level read depth across the region (plotted as a log10 aggregated over every 25 bases)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Plot conservation"),
                                                   shiny::tags$td("Tick to plot the conservation as phastCons7 across the region (aggregated over every 25 bases)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Plot constraint"),
                                                   shiny::tags$td("Tick to plot the constraint as CDTS across the region (aggregated over every 25 bases)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("SNP(s) of interest"),
                                                   shiny::tags$td("Input the coordinates of a SNP in the form chrXX:start-end of which to highlight (SNP must overlap the region to be plotted)")
                                                 )
                               ),br(),
                               h3("Plot"),
                               p("The plot is designed to help visualise a gene of interest and easily discern whether there are misannotations in potential regions of interest (e.g. overlap with a variant of unknown significance)."),
                               img(src = "ERLIN1_OMIM_reannot_example_vizER_w_help.png", 
                                   alt = "Annotated Plot",
                                   width = "100%",
                                   title = "Annotated Plot"),br(),
                               h3("ER summary table"),
                               p("The output table summarises the details of the ERs within the queried region. Each row corresponds to one ER and columns give the properties of each ER detailed below, allowing to easily query if a variant lies within any of the ERs."),
                               shiny::tags$table(id = "help-output-data",
                                                 shiny::tags$tr(
                                                   shiny::tags$th("Input name"),
                                                   shiny::tags$th("Description")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("ER_chr"),
                                                   shiny::tags$td("Chromosome expressed region is found on (1-22,  X or Y)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("ER_start"),
                                                   shiny::tags$td("Start position of expressed region (hg38)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("ER_end"),
                                                   shiny::tags$td("End position of expressed region (hg38)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("ER_width"),
                                                   shiny::tags$td("Total length in base pairs of expressed region")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Tissue"),
                                                   shiny::tags$td("The tissue from which the RNA was extracted, sequenced and analysed to derive the ER")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Mean_coverage"),
                                                   shiny::tags$td("A annotation agnostic measure of read depth (i.e. the mean number of reads overlapping each base of the ER averaged across all samples from the GTEx tissue)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Ensembl_grch38_v92_region_annot"),
                                                   shiny::tags$td("which annotation features (exon,  intron,  intergenic) the ER overlaps according to Ensembl v92")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Misannot_type"),
                                                   shiny::tags$td("The evidence by which the ER is connected to the gene of interest: 1) split_read - the ER has an overlapping split read connecting it to the gene of interest (a read with a gapped alignment to the genome), 2) overlap - the ER overlaps the gene of interest, 3) within_10Kb - the ER is within 10Kb of the gene of interest.")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Associated_gene"),
                                                   shiny::tags$td("Ensembl ID of the gene of interest")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Overlap_any_gene_v92_name"),
                                                   shiny::tags$td("Ensembl ID of the gene the expressed region overlaps")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Nearest_any_gene_v92_name"),
                                                   shiny::tags$td("Ensembl ID of the gene the ER falls closest to in terms of genomic region")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Nearest_any_gene_v92_distance"),
                                                   shiny::tags$td("The distance in bp to the nearest gene")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_annotation_type"),
                                                   shiny::tags$td("A split read intersects with the ER and either the acceptor or donor overlaps with a known exon")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_to_any_gene"),
                                                   shiny::tags$td("the Ensembl ID that the split read intersecting the ER connects to")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_ids"),
                                                   shiny::tags$td("Recount2 IDs of any split reads overlapping the ER separated by a ';'")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_count_samp"),
                                                   shiny::tags$td("start positions of any split reads overlapping the ER separated by a ';'")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_propor_samp"),
                                                   shiny::tags$td("Proportion of samples of that tissue that any split reads overlapping the ER are detected withinseparated by a ';'")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_starts"),
                                                   shiny::tags$td("start positions of any split reads overlapping the ER separated by a ';'")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Split_read_ends"),
                                                   shiny::tags$td("End positions of any split reads overlapping the ER separated by a ';'")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Mean_CDTS_percentile"),
                                                   shiny::tags$td("Mean percentile of constraint (tolerance of mutation using alignment of 7794 human genomes) across the ER (1 being the most constrained/highest evidence of functional role in humans)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Mean_phast_cons_7"),
                                                   shiny::tags$td("Mean conservation (7 species) across the ER (0 - 1 with 1 being the most conserved)")
                                                 ),
                                                 shiny::tags$tr(
                                                   shiny::tags$td("Mean_phast_cons_100"),
                                                   shiny::tags$td("Mean conservation (100 species) across the ER (0 - 1 with 1 being the most conserved)")
                                                 )
                               ),br()
                               # uiOutput("annotated")
                        ))),
             tabPanel("Contact",
                      fluidRow(
                        column(12,
                               h3("Ryten Lab"),br(), 
                               h4("This resource is generated by the Ryten Lab."),
                               p("UCL Queen Square Institute of Neurology"),
                               p("Office: 2nd floor, Russell Square House, 10-12 Russell Square, London  WC1B 5EH"),
                               a(href="https://snca.atica.um.es/", "Visit us", target="_blank"),
                               br(),br(),
                               h4("For any questions related to this resource or publication please contact:"),
                               p("Mina Ryten:",a(href="mailto:mina.ryten@ucl.ac.uk","mina.ryten@ucl.ac.uk")),
                               p("Sebastian Guefi:",a(href="mailto:manuelsebastian.guelfi@gmail.com","manuelsebastian.guelfi@gmail.com")),
                               p("David Zhang:", a(href="mailto:david.zhang.12@ucl.ac.uk","david.zhang.12@ucl.ac.uk")),
				p("Sonia Garcia-Ruiz:", a(href="mailto:s.ruiz@ucl.ac.uk","s.ruiz@ucl.ac.uk"))


                      )))
  ),
  div(class="modal")
)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ############################################################
  ##################### MAIN FUNCTION ########################
  ############################################################  
  
  genePlot <- eventReactive(input$update, {
    
    withProgress(message = 'Making plot...', value = 0.1, min = 0, max = 1, expr =  {

      ######### BODY LOADING CLASS ###########
      shinyjs::addCssClass(class = "loading", selector = "body")
      
      ######### DISABLE BUTTONS ###########
      shinyjs::disable("download_plot")
      shinyjs::disable("download_data")
      shinyjs::disable("download_bed")
      
      
      ######### STOP ZOOMING ###########
      session$sendCustomMessage(type = 'stopzoom',
                                message = list())  
      
      ######### DAVID'S FUNCTIONS ###########

      visualise_ER_example(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, 
                           txdb = ensembl_grch38_v92_genes_txdb, 
                           ensembl_gene_id_to_symbol_df = ensembl_gene_id_to_symbol_df_v92,
                           gene_id = input$geneid,
                           tissues_to_plot = input$tissue, 
                           genome_build = input$genomebuild,
                           gtex_split_read_table_mean_cov_df,
                           tissue_optimal_cut_off_max_gap_df,
                           get_constraint = input$get_constraint,
                           get_conserv = input$get_conserv,
                           get_mean_cov =  ifelse(length(input$tissue) == 1, input$get_mean_cov, FALSE),
                           propor_samples_split_read = 0.05,
                           extend_region_to_plot = ifelse(input$auto, "auto", input$extend_region_to_plot),
                           collapseTranscripts = "meta",
                           transcriptAnnotation = "gene",
                           aceview_annot = NULL,
                           add_custom_annot_track = input$add_custom_annot_track,
                           all_split_reads = F)
      
      
      data <- get_ER_table_to_display(ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific_db, 
                                      txdb = ensembl_grch38_v92_genes_txdb, 
                                      ensembl_gene_id_to_symbol_df = ensembl_gene_id_to_symbol_df_v92,
                                      gene_id = input$geneid, 
                                      tissues_to_plot = input$tissue,
                                      gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df, 
                                      extend_region_to_plot = input$extend_region_to_plot)
      ######### ENABLE BUTTONS ###########
      shinyjs::enable("download_plot")
      shinyjs::enable("download_data")
      setProgress(value = 0.99)

      ######### BODY LOADING CLASS ###########
      shinyjs::removeCssClass(class = "loading", selector = "body")
      
      ######### RETURN DATA ###########
      list(data = data)
    })
    
  })
  
  ##########################################################
  ################# VALIDATE INPUT DATA ####################
  ##########################################################
  
  validateData <- function(){
    validate(
      need(input$geneid, 'Please, type a gene.'),
      need(input$tissue, 'Please, choose a tissue.')
    )
  }
  
  
  ############################################################
  ####################### OBSERVERS ##########################
  ############################################################  
  
  observe({
	cdata <- parseQueryString(session$clientData$url_search)
	if (!is.null(cdata[['gene']])) {
		updateTextInput(session, "geneid", value = cdata[['gene']])
	}

  })


  observeEvent(input$startzoom,{
    session$sendCustomMessage(type = 'startzoom', message = list())
  })
  
  observeEvent(input$update,{
    updateTabsetPanel(session = session, inputId = "gen_browser_panel", selected = "Summary")
    updateTabsetPanel(session = session, inputId = "gen_browser_panel", selected = "Plot")
    session$sendCustomMessage(type = 'startzoom', message = list())
  })
  
  observeEvent(input$gen_browser_panel,{
    if(input$gen_browser_panel == "Plot"){
      session$sendCustomMessage(type = 'startzoom', message = list())       
    }
    else{
      session$sendCustomMessage(type = 'stopzoom', message = list())  
    }
  })
  
  observeEvent(input$auto, {
    if(input$auto){
      shinyjs::disable("extend_region_to_plot")
    }
    else
      shinyjs::enable("extend_region_to_plot")
  }, ignoreNULL = FALSE)
  
  observeEvent(input$tissue, {
    if(length(input$tissue) == 1){
      shinyjs::enable("get_mean_cov")
    }
    else{
      shinyjs::disable("get_mean_cov")
    }
  }, ignoreNULL = FALSE)
  
  
  
  
  ############################################################
  ################# BROWSER SECTION ##########################
  ############################################################
  
  ######### PLOT TAB ###########
  output$plot <- renderImage({
    validateData()
    shinyjs::disable("download_plot")
    shinyjs::disable("download_data")
    shinyjs::disable("download_bed")
    
    genePlot()
    shinyjs::enable("update")
    
    list(src = "www/OMIM_reannot_mobile_plot.png",  
         width = "95%",
         alt = "Plot",
         contentType = "image/png")
  }, deleteFile = F)
  
  ######### SUMMARY TAB ###########
  output$summary = DT::renderDataTable({
    genePlot()$data
  })
  
  ######### DOWNLOAD TAB ###########
  output$download_plot = downloadHandler(
    filename = function() {
      paste0(input$geneid,"-", input$tissue, "-", Sys.time(), "vizER_plot.png", sep="")
    },
    content = function(file) { 
      file.copy("www/OMIM_reannot_mobile_plot.png", file)
    },
    contentType = "image/png"
  )   
  output$download_data <- 
    downloadHandler(
      filename = "data.csv",
      content = function(file){
        write.csv(genePlot()$data,file)
      }
    )
  
}


onStart <- function() {
  onStop(function() {
    if(file.exists("OMIM_plot.png"))
      file.remove("OMIM_plot.png")
  })
}

# Run the application 
shinyApp(ui = ui, server = server, onStart = onStart)

