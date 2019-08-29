library(DT)
library(ggplot2)
library(tidyverse)
library(markdown)
library(shiny)
library(batman)
library(knitr)
library(shinyWidgets)

na.to.f <- function(val) { 
  return( ifelse(val == TRUE, T, F))
}

get.dr.plot <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  
  if (input$poly) {
    plt <- assay.dat %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab, group=panel_id, shape=as.factor(panel_id)))+ geom_point(size=5) + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + geom_smooth(method="lm",formula=y ~ poly(x, 5, raw=TRUE),color="red", se=F) + ggtitle('Dose-response Curve') 
    
  } else {
    plt <- assay.dat %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab, group=panel_id, shape=as.factor(panel_id)))+ geom_point(size=5) + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Dose-response Curve') 
    
  }

  return (plt)
}

get.dr.table <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  
  assay.dat <- assay.dat %>% select(auc, PAC, panel_id, prob_AIC, poly_AIC, prob_deviance) %>% unique() # conc_norm, cell_viab,

  return (assay.dat)
}

QC_filter <- function(dat) {
  dat <- dat %>% filter( !to_logical(low_PAC_flag) & !to_logical(is_within_plate_repl) & !to_logical(is_across_plate_repl) & !to_logical(across_plate_repl_flag) & !to_logical(AIC_flag) & !to_logical(DEV_flag))  # overfit_flag??? 
  return(dat)
}

get.PAC.plot <- function(input){ 
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  filt <- assay.dat %>% select(panel_id, plate_num) %>% unique()
  PAC.dat <- func.dat %>% filter( inhibitor %in% c('DMSO', 'NONE') & lab_id == input$lab_id & panel_id %in% filt$panel_id & plate_num %in% filt$plate_num) %>% mutate(plate.id = as.factor(paste(panel_id, plate_num, sep='-')))
  
  plt <- PAC.dat %>% ggplot(aes(x=cell_viab, fill=plate.id)) + geom_density(alpha=0.2) + ggtitle('Plate Controls Distribution')
  
  return(plt)
}


# ----------------------------------------------------------------------------------------------

func.dat <- read.csv('../../output/HNSCC_all_functional_data.csv', as.is=T)


server <- function(input, output, session) {

  output$dr_curve <- renderPlot({ get.dr.plot(input) })
  
  output$PAC_controls <- renderPlot({ get.PAC.plot(input)})
  
  output$DR.table <- renderDataTable({ get.dr.table(input) })
  
  output$summary <- renderPrint({ summary(func.dat) })
  
  output$table <- renderDataTable({ DT::datatable(func.dat) })
}

ui <- navbarPage("HNSCC Functional Data GUI",
                 tabPanel("Patient-Level",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput('lab_id', 'Patient ID', unique(func.dat$lab_id),
                                          selected=NULL), 
                              selectInput('inhib', 'Inhibitor', unique(func.dat$inhibitor),
                                          selected=NULL),
                              switchInput(inputId = 'poly', label = "Show poly", value = FALSE)
                            ),
                            mainPanel(
                                fluidRow(
                                  plotOutput("dr_curve"),
                                  DT::dataTableOutput("DR.table"), 
                                  plotOutput('PAC_controls')
                                  )
                              )
                          )
                 ),
                 tabPanel("Inhibitor-Level",
                          verbatimTextOutput("summary")
                 ),
                 tabPanel("About",
                          fluidRow(
                            column(12,
                                   includeMarkdown("about.md")
                            )
                          )
                 )
)


shinyApp(ui= ui, server= server)
