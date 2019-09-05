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
  
  plt <- assay.dat %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab, group=panel_id, shape=as.factor(panel_id)))+ geom_point(size=5) + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Dose-response Curve')  + theme(legend.position = "none") 
  
  if (input$poly) {
    plt <- plt  + geom_smooth(method="lm",formula=y ~ poly(x, 5, raw=TRUE),color="red", se=F) 
  }  
  
  if (input$herm) { 
    plt <- plt + geom_vline(xintercept=log10(unique(assay.dat$hermetic_transition)), color='red')  
  }
  # + geom_vline(xintercept=2)

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

get.inhib.auc.dist <- function(input) { 
  inhib.dat <- func.dat %>% filter(inhibitor == input$inhib2 & !is.na(inhibitor) ) %>% select(auc, call, lab_id, inhibitor, panel_id, plate_num) %>% unique()
  
  if (input$hist){ 
    if (input$grp) { 
      plt <- inhib.dat %>% ggplot(aes(auc, fill=call)) + geom_histogram(bins=input$bins) + ggtitle('AUC distribution')
    } else { 
      plt <- inhib.dat %>% ggplot(aes(auc)) + geom_histogram(bins=input$bins) + ggtitle('AUC distribution')
    }
  } else{
      if (input$grp) { 
        plt <- inhib.dat %>% ggplot(aes(auc, fill=call)) + geom_density(alpha=0.2) + ggtitle('AUC distribution')
      } else { 
        plt <- inhib.dat %>% ggplot(aes(auc)) + geom_density() + ggtitle('AUC distribution')
      }
  
  }
  

  return(plt)
  
}

get.atyp.plot <- function(input) { 
  inhib.dat <- func.dat %>% filter(inhibitor == input$inhib2 & !is.na(inhibitor) )
  
  plt <- inhib.dat %>% ggplot(aes(x=log10(hermetic_transition))) + geom_density(alpha=0.1, fill='red') + stat_smooth(aes(x=log10(conc_norm), y=cell_viab, group=lab_id + panel_id), geom='line', alpha=0.25, color='blue', alpha = 0.01, se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Predicted Hermetic Transitions')
  
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
  
  output$inhib_dist <- renderPlot({ get.inhib.auc.dist(input) })
  
  output$inhib_atyp <- renderPlot({ get.atyp.plot(input) })
}

ui <- navbarPage("HNSCC Functional Data GUI",
                 tabPanel("Patient-Level",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput('lab_id', 'Patient ID', unique(func.dat$lab_id),
                                          selected=NULL), 
                              selectInput('inhib', 'Inhibitor', unique(func.dat$inhibitor),
                                          selected=NULL),
                              tags$b('Display 5th-order polynomial fit'),
                              switchInput(inputId = 'poly', label = "", value = FALSE),
                              tags$b('Display predicted hermetic transition points'),
                              switchInput(inputId = 'herm', label = "", value = FALSE)
                            ),
                            mainPanel(
                                fluidRow(
                                  column(11, plotOutput("dr_curve")),
                                  column(11, DT::dataTableOutput("DR.table")), 
                                  column(11, plotOutput('PAC_controls'))
                                  )
                              )
                          )
                 ),
                 tabPanel("Inhibitor-Level",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput('inhib2', 'Inhibitor', unique(func.dat$inhibitor),
                                          selected=NULL), 
                              sliderInput("bins", "Number of Bins",
                                          min = 5, max = 30,
                                          value = 20),
                              tags$b('Inhibitor AUC distribution'),
                              switchInput(inputId = 'hist', label = "hist", value = TRUE),
                              tags$b('Group by sensitivity desgination'),
                              switchInput(inputId = 'grp', label = "", value = FALSE)
                            ),
                            mainPanel(
                              fluidRow(
                                column(11, plotOutput("inhib_dist")), 
                                column(11, plotOutput('inhib_atyp'))
                              )
                            )
                          )
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
