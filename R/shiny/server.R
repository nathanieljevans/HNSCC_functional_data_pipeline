library(DT)
library(ggplot2)
library(tidyverse)

function(input, output, session) {
  func.dat <- read.csv('../../output/HNSCC_all_functional_data.csv', as.is=T)
  
  output$plot <- renderPlot({ func.dat %>% filter(inhibitor == input$inhib & lab_id == input$lab_id) %>% ggplot(aes(x=log10(conc_norm) , y=cell_viab)) + geom_point() + geom_smooth(color='blue', se = F, method='glm', method.args=list(family=binomial(link="probit"))) + geom_smooth(color='red', se=F, method='lm', alpha=0.3) + geom_smooth(method="lm",formula=y ~ poly(x, 5, raw=TRUE),color="red", se=F) + ggtitle('Dose-response Curve') })
  
  output$summary <- renderPrint({ summary(func.dat) })
  
  output$table <- DT::renderDataTable({ DT::datatable(func.dat) })
}
