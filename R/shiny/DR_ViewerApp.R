library(DT)
library(ggplot2)
library(tidyverse)
library(markdown)
library(shiny)
library(batman)

na.to.f <- function(val) { 
  return( ifelse(val == TRUE, T, F))
}

get.dr.plot <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  
  plt <- assay.dat %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab))+ geom_point(size=5) + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + geom_smooth(method="lm",formula=y ~ poly(x, 5, raw=TRUE),color="red", se=F) + ggtitle('Dose-response Curve') 
  
  return (plt)
}

get.dr.table <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  
  assay.dat <- assay.dat %>% select(conc_norm, cell_viab, auc, PAC, prob_AIC, poly_AIC, prob_deviance) 

  return (assay.dat)
}

QC_filter <- function(dat) {
  dat <- dat %>% filter( !to_logical(low_PAC_flag) & !to_logical(is_within_plate_repl) & !to_logical(is_across_plate_repl) & !to_logical(across_plate_repl_flag) & !to_logical(AIC_flag) & !to_logical(DEV_flag))  # overfit_flag??? 
  return(dat)
}




# ----------------------------------------------------------------------------------------------

func.dat <- read.csv('../../output/HNSCC_all_functional_data.csv', as.is=T)


server <- function(input, output, session) {

  output$plot <- renderPlot({ get.dr.plot(input) })
  
  output$DR.table <- renderDataTable({ get.dr.table(input) })
  
  output$summary <- renderPrint({ summary(func.dat) })
  
  output$table <- renderDataTable({ DT::datatable(func.dat) })
}

ui <- navbarPage("HNSCC Functional Data GUI",
                 tabPanel("Plot",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput('lab_id', 'Patient ID', unique(func.dat$lab_id),
                                          selected=NULL), 
                              selectInput('inhib', 'Inhibitor', unique(func.dat$inhibitor),
                                          selected=NULL)
                            ),
                            mainPanel(
                                fluidRow(
                                  plotOutput("plot"),
                                  DT::dataTableOutput("DR.table")
                                  )
                              )
                          )
                 ),
                 tabPanel("Summary",
                          verbatimTextOutput("summary")
                 ),
                 navbarMenu("More",
                            tabPanel("Table",
                                     DT::dataTableOutput("table")
                            ),
                            tabPanel("About",
                                     fluidRow(
                                       column(6,
                                              includeMarkdown("about.md")
                                       )
                                     )
                            )
                 )
)


shinyApp(ui= ui, server= server)
