
library(MPAtools)
library(manipulate)

mpa_sim3 <- function (area, nsteps, r, pop0, K, mrate, op = FALSE, cf = NULL, closure.vec = NULL) {
    check = mpa_check(area, nsteps, r, pop0, K, mrate)
  if (check$result == 0) {
    stop(check$message)
  }
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
  if (!is.null(cf)) {
    closure = seq(cf, nsteps, by = cf)
  }
  if (!is.null(closure.vec)) {
    closure = closure.vec
  }
  if (is.null(closure.vec)) {
    closure = 0
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
    Sys.sleep(0.25)
    leaving = pop * mrate
    arriving = 0.25 * leaving[left.cells] + 0.25 * leaving[right.cells] + 
      0.25 * leaving[up.cells] + 0.25 * leaving[down.cells]
    surplus = (r * pop) * (1 - (pop/K))
    catches = u.vec * pop
    if (any(closure == i)) {
      catches = 0
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
      #png(paste(0, i, ".png", sep = ""))
      filled.contour(pop, col = rainbow(100), asp = 1, main = i, key.title = title(main = "Density\n(Org/cell)"), zlim = c(90, 810))
      #dev.off()
    }
  }
  results = list(time.series = data.frame(time = seq(0:nsteps) - 
                                            1, pop = pob/(nrows * ncols), pop.in, pop.out, catches = total.catches), 
                 pop = popi)
  return(results)
}



area <- matrix(nrow = 10, ncol = 10, 0.4)
area[4:7,3:8] <- 0

image(area)
a = mpa_sim3(area = area, nsteps = 20, r = 0.4, pop0 = 500, K = 1000, mrate = 0.1, op = T)

windows()
mpa_plot(a, "io")

###########

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

K <- 53216.2302
r <- 1.0834

manipulate(
  mpa_plot(mpa_sim2(area = area, nsteps = nsteps, pop0 = K/4, r = r, K = K, mrate = mrate, fish = fish), type = "io"),
  r=slider(0,2, initial = 1.03, step = 0.01),
  mrate=slider(0,1, initial = 0.5),
  fish = slider(0,50, initial = 10),
  nsteps = slider(1,50, initial = 50))