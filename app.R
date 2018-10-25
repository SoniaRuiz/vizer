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

source("~/R/GenomeBrowser/visualise_ER_example_gviz_v2.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  shinyjs::useShinyjs(),
  
  navbarPage("Genomic Browser",
             
             tabPanel("Browser",
                      # Application title
                      titlePanel("David Zhang et al."),
                      
                      # Sidebar with a slider input for number of bins 
                      
                      sidebarLayout(
                        sidebarPanel(
                          textInput(inputId = "geneid", label = "GeneID", value = "ENSG00000145335", width = NULL, placeholder = "Type gene ID"),
                          selectInput("tissue", "Tissue:", choices = c("Brain Cerebellar Hemisphere" = "brain_cerebellar_hemisphere",
                                                                       "Lung" = "lung",
                                                                       "Frontal Cortex" = "frontalcortexba9",
                                                                       "Adipose - subcutaneous" =	"adipose_subcutaneous",
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
                                                                       "Whole blood" =	"whole_blood"),  multiple = T, selected = "brain_cerebellar_hemisphere"),
                          numericInput(inputId = "propor_samples_split_read", label = "Propor. samples split read:", value = 0.05, step = 0.05, min = 0.05, max = 1),
                          sliderInput(inputId = "extend_region_to_plot", label = "Extend region to plot:",
                                      min = 1000, max = 50000, value = 1000, step = 500),
                          checkboxInput("get_constraint", "Get constraint", FALSE),
                          checkboxInput("get_conserv", "Get conservation", FALSE),
                          checkboxInput("all_split_reads", "All split read", FALSE),
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
                                              shiny::tags$head(shiny::tags$script(src = "jquery.elevatezoom.min.js")),
                                              singleton(
                                                shiny::tags$head(shiny::tags$script('Shiny.addCustomMessageHandler("startzoom",
                                                                                    function(message) {
                                                                                    $("#plot img").elevateZoom({
                                                                                    zoomWindowPosition:11,
                                                                                    zoomWindowWidth:550,
                                                                                    zoomWindowHeight:500,
                                                                                    scrollZoom:true
                                                                                    });
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
                            tabPanel(title = 'Summary',
                                     fluidRow(
                                       column(width = 12, 
                                              br(), br(),  
                                              DT::dataTableOutput("summary")))),
                            tabPanel(title = 'Download Data',
                                     fluidRow(
                                       column(width = 12, 
                                              br(),
                                              h2("Download page"),br(),
                                              p("Please, press next button to download the plot in PNG format:"),
                                              downloadButton(outputId = "download_plot", label = "Download Plot"), br(),br(),
                                              p("Please, press next button to download new expressed regions in CSV format:"),
                                              downloadButton(outputId = "download_data", label = "Download Data"))))
                                     )))),
             tabPanel("Paper",
                      fluidPage( 
                        titlePanel("Paper"),
                        
                        navlistPanel(
                          "Sections",
                          tabPanel("Abstract", fluidPage(
                            h1("Abstract"),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                          )),
                          tabPanel("Discussion", fluidPage(
                            h1("Discussion"),
                            p("Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi consequatur? Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam nihil molestiae consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?"),
                            p("Quod libris noster eum ne. Cetero recusabo sea at, duo laudem philosophia vituperatoribus et. Vide simul consul mea an, in audire diceret facilisi qui. An scripserit consequuntur sed, reque mandamus te duo. Pri dicunt dignissim concludaturque at, vix in eros oportere democritum, ut mel sint tation percipitur."),
                            p("Pri in suas idque, in epicurei ocurreret est. Sea ut purto omittam signiferumque, sed ex suas libris. Unum labores eloquentiam ius in. Consequat intellegat constituto vis in."),
                            p("Pro at iuvaret facilisis gubergren, sea id porro ullamcorper. Nisl meis vis cu. Mea brute fuisset te, modus deleniti et sea. Sit et autem dicit utroque, in sit justo laoreet."),
                            p("Cu sea vocibus accommodare. Ius tempor omittantur te, ex mel nisl mundi eligendi. Tempor postea animal ea vim. Ut quo sapientem mnesarchum disputando, pro an quidam patrioque.")
                          )),
                          tabPanel("Conclussions", fluidPage(
                            h1("Conclussions"),
                            p("Pri in suas idque, in epicurei ocurreret est. Sea ut purto omittam signiferumque, sed ex suas libris. Unum labores eloquentiam ius in. Consequat intellegat constituto vis in."),
                            p("Pro at iuvaret facilisis gubergren, sea id porro ullamcorper. Nisl meis vis cu. Mea brute fuisset te, modus deleniti et sea. Sit et autem dicit utroque, in sit justo laoreet."),
                            p("Cu sea vocibus accommodare. Ius tempor omittantur te, ex mel nisl mundi eligendi. Tempor postea animal ea vim. Ut quo sapientem mnesarchum disputando, pro an quidam patrioque."),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                          ))
                        ))),
             tabPanel("Help",
                      fluidRow(
                        column(12,
                               h1("README"),
                               p("This is a README detailing the columns of the data that are downloaded containing the ER definitions."),
                               p("Expressed regions (ERs) are continous segments of the genome derived from the derfinder methodolody that have evidence of being transcribed."),
                               p("ERs have been defined using RNA-sequencing data from 41 different GTEx tissues."),
                               p("ERs are connected to known genes using split reads (those with a gapped alignment to the genome) with the overarching aim of improving existing gene annotation."),
                               p("Each row represents 1 ER and each column details one property of the corresponding ER."),br(),
                               h1("Annotated Properties"),
                               tableOutput('annotated_properties'),br(),
                               h1("Annotated Plot"),
                               uiOutput("annotated")
                        ))
             ))
  )

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ############################################################
  ##################### MAIN FUNCTION ########################
  ############################################################  
  
  gene_plot <- eventReactive(input$update, {
    
    withProgress(message = 'Making plot...', value = 0.1, min = 0, max = 1, expr =  {
      
      ######### DISABLE BUTTONS ###########
      shinyjs::disable("update")
      shinyjs::disable("download_plot")
      shinyjs::disable("download_data")
      
      ######### STOP ZOOMING ###########
      session$sendCustomMessage(type = 'stopzoom',
                                message = list())  
      
      ######### DAVID'S FUNCTIONS ###########
      visualise_ER_example(ERs_w_annotation_df = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, 
                           txdb = ensembl_grch38_v92_genes_txdb, 
                           ensembl_gene_id_to_symbol_df = ensembl_gene_id_to_symbol_df_v92,
                           gene_id = input$geneid,
                           tissues_to_plot = input$tissue, 
                           genome_build = input$genomebuild,
                           gtex_split_read_table_mean_cov_df,
                           tissue_optimal_cut_off_max_gap_df,
                           get_constraint = input$get_constraint,
                           get_conserv = input$get_conserv,
                           get_mean_cov = F,
                           propor_samples_split_read = input$propor_samples_split_read,
                           extend_region_to_plot = input$extend_region_to_plot,
                           collapseTranscripts = "meta",
                           transcriptAnnotation = "gene",
                           aceview_annot = NULL,
                           add_custom_annot_track = NULL,
                           all_split_reads = input$all_split_reads)
      dev.print(file = "OMIM_plot.png", device = png, res = 600, width = 10, height = 11.69/2, units = "in")
      data <- get_ER_table_to_display(ERs_w_annotation_df = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, 
                                      txdb = ensembl_grch38_v92_genes_txdb, 
                                      ensembl_gene_id_to_symbol_df = ensembl_gene_id_to_symbol_df_v92,
                                      gene_id = input$geneid, 
                                      tissues_to_plot = input$tissue,
                                      gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df, 
                                      extend_region_to_plot = input$extend_region_to_plot)
      
      ######### ENABLE BUTTONS ###########
      shinyjs::enable("update")
      shinyjs::enable("download_plot")
      shinyjs::enable("download_data")
      
      ######### RETURN DATA ###########
      list(data = data)
    })
  })
  
  
  ############################################################
  ####################### OBSERVERS ##########################
  ############################################################  
  
  observeEvent(input$startzoom,{
    session$sendCustomMessage(type = 'startzoom', message = list())
  })
  
  observeEvent(input$update,{
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
  
  ############################################################
  ################# BROWSER SECTION ##########################
  ############################################################
  
  ######### PLOT TAB ###########
  output$plot <- renderImage({
    shinyjs::disable("download_plot")
    shinyjs::disable("download_data")
    gene_plot()
    list(src = "OMIM_plot.png",  
         width = "95%",
         alt = "Plot",
         contentType = "image/png")
  }, deleteFile = FALSE)
  
  ######### SUMMARY TAB ###########
  output$summary = DT::renderDataTable({
    gene_plot()$data
  })
  
  ######### DOWNLOAD TAB ###########
  output$download_plot = downloadHandler(
    filename = function() {
      paste0("plot-", Sys.time(), ".png", sep="")
    },
    content = function(file) {
      dev.copy(file = file, device = png, res = 600, width = 10, height = 11.69/2, units = "in")
      dev.off()
    })   
  output$download_data <- 
    downloadHandler(
      filename = "data.csv",
      content = function(file){
        write.csv(gene_plot()$data,file)
      }
    )
  
  ############################################################
  ################# README SECTION ###########################
  ############################################################
  
  output$annotated_properties <- renderTable({
    tb <- read.csv(file = "./www/ERproperties.csv", header=T, sep=",")
    return(tb)
  })
  output$annotated <- renderUI({
    img(src = "annotated_plot.png", 
        width = "100%",
        alt = "Annotated Plot")
  })
}
# Run the application 
shinyApp(ui = ui, server = server)