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

source("~/R/GenomeBrowser/visualise_ER_example_gviz_sonia.R")

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
                                     fluidRow(
                                       column(width = 12, 
                                              br(), br(),  
                                              plotOutput("plot")
                                       )
                                     )
                            ),
                            tabPanel(title = 'Summary',
                                     fluidRow(
                                       column(width = 12, 
                                              br(), br(),  
                                              DT::dataTableOutput("summary")
                                       )
                                     )
                            ),
                            tabPanel(title = 'Download Data',
                                     fluidRow(
                                       column(width = 12, 
                                              br(),
                                              h2("Download page"),br(),
                                              p("Please, press next button to download the plot in PNG format:"),
                                              downloadButton(outputId = "download_plot", label = "Download Plot"), br(),br(),
                                              p("Please, press next button to download new expressed regions in CSV format:"),
                                              downloadButton(outputId = "download_data", label = "Download Data")
                                       )
                                     )
                            )#,
                            #tabPanel("Help", htmlOutput("help"))
                          )
                        )
                      )
                    ),
             tabPanel("Paper",
                      titlePanel("David Zhang et al."),
                      sidebarLayout(
                        sidebarPanel(
                          h3("Details"),
                          h4(code("Title:"), " ..."), 
                          h4(code("Author:"), " David Zhang"), 
                          h4(code("Publication:"), " ..."),
                          h4(code("Volume:"), " ..."),
                          h4(code("Issue:"), " ..."),
                          h4(code("Year:"), " 2018")
                        ),
                        mainPanel(
                          h1("David Zhang et. al"),
                          h3("Abstract"),
                          p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."),
                          h3("Discussion"),
                          p("Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi consequatur? Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam nihil molestiae consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?")
                        )
                      )),
             tabPanel("Help",
                      fluidPage(
                        titlePanel("Help"),
                        
                        navlistPanel(
                          "How to use",
                          tabPanel("Get constraint", fluidPage(
                            h1("Get constraint"),
                            h3("What does it means?"),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                          )),
                          tabPanel("Get conservation", fluidPage(
                            h1("Get conservation"),
                            h3("What does it means?"),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                          )),
                          "Results",
                          tabPanel("New expressed regions", fluidPage(
                            h1("New expressed regions"),
                            h3("Intron"),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."),
                            h3("Intergenic"),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                            
                          )),
                          tabPanel("Interpreting the result", fluidPage(
                            h1("What does it means?"),
                            p("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
                          ))
                        )
                      )
  ))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  observeEvent(input$update,{
    updateTabsetPanel(session = session, inputId = "gen_browser_panel", selected = "Plot")})
  

  
  gene_plot <- eventReactive(input$update, {
   
    withProgress(message = 'Making plot...', value = 0.1, min = 0, max = 1, expr =  {
      
      shinyjs::disable("update")
      shinyjs::disable("download_plot")
      shinyjs::disable("download_data")
      
      visualise_ER_example(ERs_w_annotation_df = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, 
                           txdb = ensembl_grch38_v92_genes_txdb, 
                           ensembl_gene_id = input$geneid,
                           tissues_to_plot = input$tissue, 
                           genome_build = input$genomebuild,
                           gtex_split_read_table_mean_cov_df,
                           tissue_optimal_cut_off_max_gap_df,
                           get_constraint = input$get_constraint,
                           get_conserv = input$get_conserv,
                           propor_samples_split_read = input$propor_samples_split_read,
                           extend_region_to_plot = input$extend_region_to_plot,
                           collapseTranscripts = "meta",
                           transcriptAnnotation = "gene",
                           aceview_annot = NULL,
                           add_custom_annot_track = NULL,
                           all_split_reads = input$all_split_reads)
      
      data <- get_ER_table_to_display(ERs_w_annotation_df = ERs_w_annotation_all_tissues_width_ab_3_no_cells_sex_specific, 
                                      txdb = ensembl_grch38_v92_genes_txdb, 
                                      ensembl_gene_id_to_symbol_df = input$geneid,
                                      gene_id = input$geneid, 
                                      tissues_to_plot = input$tissue,
                                      gtex_split_read_table_mean_cov_df = gtex_split_read_table_mean_cov_df, 
                                      extend_region_to_plot = input$extend_region_to_plot)
      shinyjs::enable("update")
      shinyjs::enable("download_plot")
      shinyjs::enable("download_data")
      
      list(data = data)
    })
  })
  

  output$plot <- renderPlot({
    shinyjs::disable("download_plot")
    shinyjs::disable("download_data")
    gene_plot()
  })
  output$summary = DT::renderDataTable({
    gene_plot()$data
  })
  outputOptions(output, "plot", suspendWhenHidden = FALSE)
  outputOptions(output, "summary", suspendWhenHidden = FALSE)
  
  output$download_plot = downloadHandler(
    filename = function() {
      paste0("plot-", Sys.Date(), ".png", sep="")
    },
    content = function(file) {
      dev.copy(file = file, device = png, res = 600, width = 10, height = 11.69/2, units = "in")
      dev.off()
    })   
  output$download_data <- 
    downloadHandler(
      filename = "data.csv",
      content = function(file){
        write.csv(head(gene_plot()$data, n=6),file)
      }
    )
  #output$help <- renderText({
  #  HTML("<br/><h2>Help to the user</h2><br/><p>We can add here some useful information to the user.</p><br/>")
  #})
  #outputOptions(output, "help", suspendWhenHidden = FALSE)
}
# Run the application 
shinyApp(ui = ui, server = server)