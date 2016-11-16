#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   
   # Application title
   titlePanel("Optimum MPA size and effort"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("u",
                     "Proportion of organisms:",
                     min = 0,
                     max = 1,
                     value = 0.5),
         sliderInput("nsteps",
                     "Number of iterations",
                     min = 0, max = 100,
                     value = 50),
         sliderInput("r",
                     "Intrinsic populatin growth rate",
                     min = 0,
                     max = 2, 
                     value = 1.0834)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
))

# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
  #build area here
  
  rev <- function(area, nsteps, pop0, r, K, mrate, closure.vec = NULL){
    sim1 <- mpa_sim(area = area, nsteps = nsteps, pop0 = K/2, r = r, K = K, mrate = mrate, closure.vec = closure.vec, op = F)
    t <- sim1$time.series$time
    c <- sim1$time.series$catches
    
    rev = sum(((c*20)-730)/((1+0.05)^t))
    
    return(rev)
  }
  
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      hist(x, breaks = bins, col = 'darkgray', border = 'white')
   })
})

# Run the application 
shinyApp(ui = ui, server = server)

