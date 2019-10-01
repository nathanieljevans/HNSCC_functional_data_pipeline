library(DT)
library(ggplot2)
library(tidyverse)
library(markdown)
library(shiny)
library(batman)
library(knitr)
library(shinyWidgets)


get.dr.plot <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  
  print(assay.dat)
  
  plt <- assay.dat %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab, group=panel_id, shape=as.factor(panel_id)))+ geom_point(size=5) + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Dose-response Curve')  + theme(legend.position = "none") 
  
  if (input$poly) {
    plt <- plt  + geom_smooth(method="lm",formula=y ~ poly(x, 5, raw=TRUE),color="red", se=F) 
  }  
  
  if (input$herm) { 
    assay.dat %>% select(conc_norm, atyp_prob) %>% print(.)
    plt <- plt + geom_col(aes(x=log10(conc_norm), y=as.numeric(atyp_prob), color='red'),alpha=0.3)
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
  
  if (input$grp) { 
    plt <- inhib.dat %>% ggplot(aes(x=log10(hermetic_transition))) + geom_density(alpha=0.05, fill='red') + stat_smooth(aes(x=log10(conc_norm), y=cell_viab, group=lab_id + panel_id, color=call), geom='line', alpha=0.5, alpha = 0.01, se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Predicted Hermetic Transitions')
  } else {
    plt <- inhib.dat %>% ggplot(aes(x=log10(hermetic_transition))) + geom_density(alpha=0.05, fill='red') + stat_smooth(aes(x=log10(conc_norm), y=cell_viab, group=lab_id + panel_id), geom='line', alpha=0.25, color='blue', alpha = 0.01, se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Predicted Hermetic Transitions')
  }
  return(plt)
}


get.pat.sens.plot <- function(input) { 
  
  assay.dat <- func.dat %>% QC_filter(.) %>% select(lab_id, inhibitor, auc, call) %>% unique() %>% group_by(lab_id, inhibitor) %>% summarize(auc = mean(auc)) %>% ungroup() %>% data.frame()
  
  assay.dat <- assay.dat %>% group_by(inhibitor) %>% mutate(AUC_z.score = scale(auc)) %>% ungroup()
  
  if (input$sens){
    dat <- assay.dat %>% filter(lab_id == input$pat2) %>% arrange(AUC_z.score) %>% head(input$naucs)
  } else{
    dat <- assay.dat %>% filter(lab_id == input$pat2) %>% arrange(desc(AUC_z.score)) %>% head(input$naucs)
  }
  inhib_order <- dat$inhibitor
  
  plt <- dat %>% ggplot(aes(x=factor(inhibitor, level=inhib_order), y=AUC_z.score)) + geom_col(alpha=0.5) + theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle('Sample Inhibitor Sensitivity')
  
  return(plt)
}

get.pat.sens.tab <- function(input) { 
  assay.dat <- func.dat %>% filter(lab_id == input$pat2) %>% QC_filter(.) %>% select(lab_id, inhibitor, plate_num, panel_id, auc, call) %>% unique()#%>% group_by(lab_id, inhibitor) %>% summarize(auc = mean(auc)) %>% ungroup() %>% data.frame()
  #print(assay.dat)
  if (input$sens) {
    tab <- assay.dat %>% filter(call == 'sens') %>% select(lab_id, panel_id, inhibitor, auc, call) %>% arrange(auc)
  } else {
    tab <- assay.dat %>% filter(call == 'res') %>% select(lab_id, panel_id, inhibitor, auc, call) %>% arrange(desc(auc))
  }
  return(tab)
}


# ----------------------------------------------------------------------------------------------

func.dat <- read.csv('../../output/HNSCC_all_functional_data.csv', as.is=T) %>% mutate(atyp_prob = as.numeric(atyp_prob))

server <- function(input, output, session) {

  output$dr_curve <- renderPlot({ get.dr.plot(input) })
  
  output$PAC_controls <- renderPlot({ get.PAC.plot(input)})
  
  output$DR.table <- renderDataTable({ get.dr.table(input) })
  
  output$summary <- renderPrint({ summary(func.dat) })
  
  output$table <- renderDataTable({ DT::datatable(func.dat) })
  
  output$inhib_dist <- renderPlot({ get.inhib.auc.dist(input) })
  
  output$inhib_atyp <- renderPlot({ get.atyp.plot(input) })
  
  output$pat_sens_plot <- renderPlot({ get.pat.sens.plot(input) })
  
  output$pat_sens_tab <- renderDataTable({ get.pat.sens.tab(input) })
}

ui <- navbarPage("HNSCC Functional Data GUI",
                 tabPanel("Assay-Level",
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
                 tabPanel("Patient-Level",
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("naucs", "Number of inhibitors to display",
                                                     min = 5, max = 100,
                                                     value = 30),
                              selectInput('pat2', 'Inhibitor', unique(func.dat$lab_id),
                                          selected=NULL), 
                              tags$b('Sort by resistant or sensitive'),
                              switchInput(inputId = 'sens', label = "sens", value = TRUE)
                            ),
                            mainPanel(
                              fluidRow(
                                column(11, plotOutput("pat_sens_plot")), 
                                column(11, DT::dataTableOutput("pat_sens_tab"))
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


shinyApp(ui=ui, server=server)
