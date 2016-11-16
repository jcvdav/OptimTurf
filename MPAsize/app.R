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
                     step = 0.1,
                     value = 0.5),
         sliderInput("nsteps",
                     "Number of iterations",
                     min = 0,
                     max = 100,
                     value = 50),
         sliderInput("r",
                     "Intrinsic populatin growth rate",
                     min = 0,
                     max = 2,
                     step = 0.1,
                     value = 1.0834),
         sliderInput("K",
                     "Carrying capacity",
                     min = 0,
                     max = 1,
                     step = 0.01,
                     value = 0.053),
         sliderInput("mrate",
                     "Movement rate",
                     min = 0,
                     max = 1,
                     step = 0.1,
                     value = 0.5)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
))

# Define server logic
server <- shinyServer(function(input, output) {
  
  rev <- function(area, nsteps, pop0, r, K, mrate, closure.vec = NULL){
    sim1 <- mpa_sim(area = area, nsteps = nsteps, pop0 = K/2, r = r, K = K, mrate = mrate, closure.vec = closure.vec, op = F)
    t <- sim1$time.series$time
    c <- sim1$time.series$catches
    
    rev = sum(((c*20)-730)/((1+0.05)^t))
    
    return(rev)
  }
  
  #build area here
  
  u <- reactive({input$u})
  r <- reactive({input$r})
  K <- reactive({input$K})
  mrate <- reactive({input$mrate})
  nsteps <- reactive({input$nsteps})
  
  area <- matrix(nrow = 10, ncol = 10, 0.5)
  area[3:6,3:7] <- 0
  
  u.vec = area
    pop <- matrix(nrow = 10, ncol = 10, K/2)
    popi = array(NA, dim = c(10, 10, 50))
    left.cells = c(ncols, 1:(ncols - 1))
    right.cells = c(2:ncols, 1)
    up.cells = c(nrows, 1:(nrows - 1))
    down.cells = c(2:nrows, 1)
    
    inside = u.vec == 0
    outside = u.vec > 0
    ones = matrix(1, nrow = nrows, ncol = ncols)
    pop.in = sum(pop[inside], na.rm = TRUE)/sum(ones[inside], 
        na.rm = TRUE)
    pop.out = sum(pop[outside], na.rm = TRUE)/sum(ones[outside], 
        na.rm = TRUE)
    pob = sum(pop, na.rm = TRUE)
    total.catches = sum(u.vec * pop, na.rm = TRUE)
    if (op == TRUE) {
        windows()
        filled.contour(pop, color = terrain.colors, asp = 1, 
            main = "0", key.title = title(main = "Density\n(Org/cell)"))
    }
  
  output$distPlot <- renderPlot({

    for (i in 1:nsteps) {
        leaving = pop * mrate
        arriving = 0.25 * leaving[left.cells] + 0.25 * leaving[right.cells] + 0.25 * leaving[up.cells] + 0.25 * leaving[down.cells]
        surplus = (r * pop) * (1 - (pop/K))
        catches = u.vec * pop
        pop = pop + surplus - catches - leaving + arriving
        popi[, , i] = pop
        pop.in[i + 1] = sum(pop[inside], na.rm = TRUE)/sum(ones[inside], 
            na.rm = TRUE)
        pop.out[i + 1] = sum(pop[outside], na.rm = TRUE)/sum(ones[outside], 
            na.rm = TRUE)
        pob[i + 1] = sum(pop, na.rm = TRUE)
        total.catches[i + 1] = sum(catches, na.rm = TRUE)
        if (op == TRUE) {
            filled.contour(pop, col = rainbow(100), asp = 1, 
                main = i, key.title = title(main = "Density\n(Org/cell)"))
        }
    }
    results = list(time.series = data.frame(time = seq(0:nsteps) - 
        1, pop = pob/(nrows * ncols), pop.in, pop.out, catches = total.catches), 
        pop = popi)
    return(results)
  })
   
   # output$distPlot <- renderPlot({
   #    # generate bins based on input$bins from ui.R
   #    x    <- faithful[, 2] 
   #    bins <- seq(min(x), max(x), length.out = input$bins + 1)
   #    
   #    # draw the histogram with the specified number of bins
   #    hist(x, breaks = bins, col = 'darkgray', border = 'white')
   # })
})

# Run the application 
shinyApp(ui = ui, server = server)

