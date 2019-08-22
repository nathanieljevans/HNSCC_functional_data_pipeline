library(markdown)
library(shiny)

func.dat <- read.csv('../../output/HNSCC_all_functional_data.csv', as.is=T)

navbarPage("HNSCC Functional Data GUI",
           tabPanel("Plot",
                    sidebarLayout(
                      sidebarPanel(
                        selectInput('lab_id', 'Patient ID', unique(func.dat$lab_id),
                                    selected=NULL), 
                        selectInput('inhib', 'Inhibitor', unique(func.dat$inhibitor),
                                    selected=NULL)
                        ),
                      mainPanel(plotOutput("plot"))
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
