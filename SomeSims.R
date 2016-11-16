## Estimate r and K parameters for Lobster in Natividad

library(tidyverse)
library(MPAtools)

rev <- function(area, nsteps, pop0, r, K, mrate, closure.vec = NULL){
  sim1 <- mpa_sim(area = area, nsteps = nsteps, pop0 = K/2, r = r, K = K, mrate = mrate, closure.vec = closure.vec, op = F)
  t <- sim1$time.series$time
  c <- sim1$time.series$catches
  
  rev = sum(((c*20)-730)/((1+0.05)^t))
  
  return(rev)
}


mpa_sim2 <- function (area, nsteps, r, pop0, K, mrate, op = FALSE, cf = NULL, fish=1) {
  u.vec = area
  size = dim(u.vec)
  nrows = size[1]
  ncols = size[2]
  pop <- matrix(nrow = nrows, ncol = ncols)
  popi = array(NA, dim = c(nrows, ncols, nsteps))
  pop[] = pop0
  left.cells = c(ncols, 1:(ncols - 1))
  right.cells = c(2:ncols, 1)
  up.cells = c(nrows, 1:(nrows - 1))
  down.cells = c(2:nrows, 1)
  
  if (fish > 0){
  fish = seq(fish, nsteps, by = fish)
  }
  
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
  for (i in 1:nsteps) {
    leaving = pop * mrate
    arriving = 0.25 * leaving[left.cells] + 0.25 * leaving[right.cells] + 
      0.25 * leaving[up.cells] + 0.25 * leaving[down.cells]
    surplus = (r * pop) * (1 - (pop/K))
    catches = u.vec * pop
    if (any(fish == i)) {
      u.vec2 = matrix(nrow = 10, ncol = 10, u.vec[1])
      catches = u.vec2 * pop
    }
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
}

rev2 <- function(area, nsteps, pop0, r, K, mrate, fish){
  sim1 <- mpa_sim2(area = area, nsteps = nsteps, pop0 = K/2, r = r, K = K, mrate = mrate, fish = fish)
  t <- sim1$time.series$time
  c <- sim1$time.series$catches
  
  rev = sum(((c*20)-730)/((1+0.05)^t))
  
  return(rev)
}

###

## Optimum size
K <- 53216.2302
r <- 1.0834

area <- matrix(nrow =10, ncol = 10, 0.5)
area[5,5] = 0 

rev(area = area, nsteps = 50, pop0 = K/2, r = r, K = K, mrate = 0.1)

size <- seq(1:100)
u <- seq(0.1, 1, by = 0.1)

Rev <- matrix(nrow = length(u), ncol = length(size), 0)

for (i in 1:length(u)){
  for (j in 1:100){
    area <- matrix(nrow = 10, ncol = 10, u[i])
    area[1:j] <- 0
    
    Rev[i,j] <- rev(area = area, nsteps = 50, pop0 = K/2, r = r, K = K, mrate = 0.1)
    
  }
}

Catch_Rate <- u
Percentage_of_Closure<- size

image(Catch_Rate, Percentage_of_Closure, Rev)
contour.default(u, size, Rev, add = T)


#Optimum closure times

Rev2 <- matrix(nrow = 51, ncol = 1)
Time <- seq(0,50)

area <- matrix(nrow = 10, ncol = 10, 0.6)
area[4:6, 4:6] <- 0
area[6,7] <- 0

for (j in 0:50){
  
  Rev2[j+1] <- rev2(area = area, nsteps = 50, pop0 = K/2, r = r, K = K, mrate = 0.1, fish = j)

}

plot(Time, Rev2)


manipulate(
  mpa_plot(mpa_sim2(area = area, nsteps = nsteps, pop0 = K/2, r = r, K = K, mrate = mrate, fish = fish), type = "io"),
  r=slider(0,2, initial = 1.03, step = 0.01),
  mrate=slider(0,1, initial = 0.5),
  fish = slider(0,50, initial = 10),
  nsteps = slider(1,50, initial = 50))

mpa_plot(mpa_sim2(area = area, nsteps = 50, pop0 = K/2, r = r, K = K, mrate = 0.5, fish = 1), type = "io")

mpa_plot(mpa_sim2(area = area, nsteps = 50, pop0 = K/2, r = r, K = K, mrate = 0.5, fish = 1), type = "c")


