library(DT)
library(ggplot2)
library(tidyverse)
library(markdown)
library(shiny)
library(batman)
library(knitr)
library(shinyWidgets)
library(shinythemes)

# -------------------------------------------------------------------------------------------------------------------
######## GLOBALS ####################################################################################################
# -------------------------------------------------------------------------------------------------------------------
OUR.THEME = "flatly"           # shiny theme 
POLY.FIT.ORDER = 5             # order of "overfitted" regression on Page 1 Dose-Response Plot
# -------------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------------
# Page 1: Top plot, Dose-Response plots, stratified by sensitivity 
# -------------------------------------------------------------------------------------------------------------------
get.dr.plot <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  print(head(assay.dat))
  
  plt <- assay.dat %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab, group=panel_id, shape=as.factor(panel_id)))+ geom_point(size=5) + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Dose-response Curve')  + theme(legend.position = "none") 
  
  if (input$poly) {
    plt <- plt  + geom_smooth(method="lm",formula=y ~ poly(x, POLY.FIT.ORDER, raw=TRUE),color="red", se=F) 
  }  
  
  if (input$herm) { 
    assay.dat %>% select(conc_norm, atyp_prob) %>% print(.)
    plt <- plt + geom_col(aes(x=log10(conc_norm), y=atyp_prob, color='red'),alpha=0.02)
  }

  return (plt)
}


# -------------------------------------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------------------------------------
get.dr.table <- function(input){ 
  
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  
  assay.dat <- assay.dat %>% select(auc, PAC, panel_id, prob_AIC, poly_AIC, prob_deviance) %>% unique() # conc_norm, cell_viab,

  return (assay.dat)
}


# -------------------------------------------------------------------------------------------------------------------
# QUALITY CONTROL FUNCTION
# -------------------------------------------------------------------------------------------------------------------
# This function is paramount in removing spurrious or duplicate observations so as not to confound visualization. 
QC_filter <- function(dat) {
  true = c('TRUE')
  false = c('FALSE','NA', '', NA)
  dat <- dat %>% filter( low_PAC_flag %in% false & 
                           is_within_plate_repl %in% false &
                           is_across_plate_repl %in% false & 
                           across_plate_repl_flag %in% false  & 
                           AIC_flag %in% false & 
                           DEV_flag %in% false)  
 
  return(dat)
}


# -------------------------------------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------------------------------------
get.PAC.plot <- function(input){ 
  assay.dat <- func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% QC_filter(.) 
  filt <- assay.dat %>% select(panel_id, plate_num) %>% unique()
  PAC.dat <- func.dat %>% filter( inhibitor %in% c('DMSO', 'NONE') & lab_id == input$lab_id & panel_id %in% filt$panel_id & plate_num %in% filt$plate_num) %>% mutate(plate.id = as.factor(paste(panel_id, plate_num, sep='-')))
  
  plt <- PAC.dat %>% ggplot(aes(x=cell_viab, fill=plate.id)) + geom_density(alpha=0.2) + ggtitle('Plate Controls Distribution')
  
  return(plt)
}


# -------------------------------------------------------------------------------------------------------------------
# PAGE 3: AUC DISTRIBUTION 
# -------------------------------------------------------------------------------------------------------------------
# Generates the auc distribution, and stratfies by sensitivity; Optional to use histogram/density by user choice.
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


# -------------------------------------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------------------------------------
get.atyp.plot <- function(input) { 
  inhib.dat <- func.dat %>% filter(inhibitor == input$inhib2 & !is.na(inhibitor) )
  
  if (input$grp) { 
    plt <- inhib.dat %>% ggplot(aes(x=log10(conc_norm), y=atyp_prob, group=log10(conc_norm))) + geom_violin(alpha=0.05, fill='red') + stat_smooth(aes(x=log10(conc_norm), y=cell_viab, group=lab_id + panel_id, color=call), geom='line', alpha = 0.01, se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Predicted Atypical Transitions')
  } else {
    plt <- inhib.dat %>% ggplot(aes(x=log10(conc_norm), y=atyp_prob, group=log10(conc_norm))) + geom_violin(alpha=0.05, fill='red') + stat_smooth(aes(x=log10(conc_norm), y=cell_viab, group=lab_id + panel_id), geom='line', color='blue', alpha = 0.01, se = F, method='glm', method.args=list(family=binomial(link="probit"))) + ggtitle('Predicted Atypical Transitions')
  }
  return(plt)
}



# -------------------------------------------------------------------------------------------------------------------
# Page 2: Z-score Ranked Inhibitors
# -------------------------------------------------------------------------------------------------------------------
# 
get.pat.sens.plot <- function(input) { 
  
  assay.dat <- func.dat %>% QC_filter(.) %>% select(lab_id, inhibitor, auc, call) %>%
          unique() %>% group_by(lab_id, inhibitor) %>% summarize(auc = mean(auc)) %>% 
          ungroup() %>% data.frame() %>% group_by(inhibitor) %>% 
          mutate(AUC_z.score = scale(auc)) %>% ungroup()
  
  this.theme <- theme(
          axis.text.x = element_text(color = "black", size = 12, angle = 90, hjust = 1, vjust = 1, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))
  
  if (input$sens){
    dat <- assay.dat %>% filter(lab_id == input$pat2) %>% arrange(AUC_z.score) %>% head(input$naucs) 
  } else{
    dat <- assay.dat %>% filter(lab_id == input$pat2) %>% arrange(desc(AUC_z.score)) %>% head(input$naucs)
  }
  inhib_order <- dat$inhibitor
  
  plt <- dat %>% ggplot(aes(x=factor(inhibitor, level=inhib_order), y=AUC_z.score)) + geom_col(alpha=0.5) + 
                  this.theme + ggtitle('Sample Inhibitor Sensitivity') + xlab('Drug Name') + ylab('Relative Sensitivity (AUC Z-score)')
  
  return(plt)
}


# -------------------------------------------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------------------------------------------
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


# -------------------------------------------------------------------------------------------------------------------
################################## DATA QUALITY CONTROL AND PREPROCESSING ###########################################
# -------------------------------------------------------------------------------------------------------------------
func.dat <- read.csv('../../output/HNSCC_all_functional_data.csv', as.is=T) %>% 
                                  mutate(atyp_prob = as.numeric(atyp_prob)) %>% QC_filter()


# -------------------------------------------------------------------------------------------------------------------
################################## [SERVER]  OUTPT DESIGNATION  [SERVER] ############################################
# -------------------------------------------------------------------------------------------------------------------
server <- function(input, output, session) {

  
                                    # ---------------------------------------------------
                                    ############## PAGE 1: ASSAY LEVEL ##################
                                    # ---------------------------------------------------
  # Dose response curve; stratified by sensitivity label
  output$dr_curve <- renderPlot({ get.dr.plot(input) })
  
  # 
  output$PAC_controls <- renderPlot({ get.PAC.plot(input)})
  
  output$DR.table <- renderDataTable({ get.dr.table(input) })
  
  
                                    # ---------------------------------------------------
                                    ############## PAGE 2: PATIENT LEVEL ################
                                    # ---------------------------------------------------
  output$summary <- renderPrint({ summary(func.dat) })
  
  output$table <- renderDataTable({ DT::datatable(func.dat) })
  
  
                                    # ---------------------------------------------------
                                    ############## PAGE 3: INHIBITOR LEVEL ##############
                                    # ---------------------------------------------------
  output$inhib_dist <- renderPlot({ get.inhib.auc.dist(input) })
  
  output$inhib_atyp <- renderPlot({ get.atyp.plot(input) })
  
  
                                    # ---------------------------------------------------
                                    ############## PAGE 4: ABOUT LEVEL ##################
                                    # ---------------------------------------------------
  # Z-score ranked inhibitor plot 
  output$pat_sens_plot <- renderPlot({ get.pat.sens.plot(input) }, height=500)
  
  output$pat_sens_tab <- renderDataTable({ get.pat.sens.tab(input) })
}


# -------------------------------------------------------------------------------------------------------------------
#####################################################################################################################
###################################### APP STRUCTURE AND DESIGN #####################################################
#####################################################################################################################
# -------------------------------------------------------------------------------------------------------------------
ui <- navbarPage("HNSCC Functional Data GUI",
                 
                 # ------------------------------------------------------------------------------------------------------------
                 ############################################# ASSAY LEVEL PAGE ###############################################
                 # ------------------------------------------------------------------------------------------------------------
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
                 
                 # ------------------------------------------------------------------------------------------------------------
                 ################################################ Patient Level Page  #########################################
                 # ------------------------------------------------------------------------------------------------------------
                 tabPanel("Patient-Level",
                          sidebarLayout(
                            sidebarPanel(
                              sliderInput("naucs", "Number of inhibitors to display",
                                                     min = 5, max = 50,
                                                     value = 20),
                              selectInput('pat2', 'Inhibitor', unique(func.dat$lab_id),
                                          selected=NULL), 
                              tags$b('Sort by resistant or sensitive'),
                              switchInput(inputId = 'sens', label = "sens", value = TRUE)
                            ),
                            mainPanel(
                              fluidRow(
                                column(width = 11, offset = 0, style='padding-left:0px; padding-right:0px; 
                                       padding-top:5px; padding-bottom:150px', plotOutput("pat_sens_plot")), 
                                column(11, DT::dataTableOutput("pat_sens_tab"))
                              )
                            )
                          )
                 ),
                 
                 # ---------------------------------------------------------------------------------------------------------
                 ################################################ INHIBITOR PAGE ###########################################
                 # ---------------------------------------------------------------------------------------------------------
                 tabPanel("Inhibitor-Level",
                          sidebarLayout(
                            sidebarPanel(
                              selectInput('inhib2', 'Inhibitor', unique(func.dat$inhibitor),
                                          selected=NULL), 
                              sliderInput("bins", "Number of Bins",
                                          min = 5, max = 30,
                                          value = 20),
                              tags$b('Histogram/Density'),
                              switchInput(inputId = 'hist', label = "", value = TRUE),
                              tags$b('Include Sensitivity'),
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
                 
                 # -------------------------------------------------------------------------------------------------------------
                 ################################################ ABOUT PAGE ###################################################
                 # -------------------------------------------------------------------------------------------------------------
                 tabPanel("About",
                          fluidRow(
                            column(12,
                                   includeMarkdown("about.md")
                            )
                          )
                 )
)


# -------------------------------------------------------------------------------------------------------------------
################################################ APP CALL ###########################################################
# -------------------------------------------------------------------------------------------------------------------
shinyApp(ui=fluidPage(theme = shinytheme(OUR.THEME), ui), 
         server=server)
