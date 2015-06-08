########################################################
# The MIT License (MIT)
#
# Copyright (c) 2014 Florian D. Schneider & Sonia Kéfi
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
########################################################

######################
## mapping function ##
######################

# required to vectorise the counting and plotting. 
# returns a map of the landscape to translate it into a vector with boundaries and another one to back-translate it to a vector without boundaries into the global environment. Needs to be called only once for the dimensions of the lattice. 


mapping <- function(width, height, boundary = "periodic", i_matrix = matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)) {
  
  # derive helper vectors for counting: 
  # transformation vector for evaluation at the border of the grid
  # set evaluation matrix 
  X <- matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
  # setting the border of the evaluation matrix X
  X <- cbind(X[,width], X, X[,1] )  
  X <- rbind(X[height,], X, X[1,] ) 
  # transformation vector which adds the border to the lattice:
  x_with_border <- as.integer(t(X))
  
  assign("x_with_border", as.integer(t(X))  , envir = .GlobalEnv )
  
  # from the matrix X (with copied border cells), which cells are the actual cells (reverse transformation vector of the previous lines) 
  #x_to_evaluate <- sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )    
  
  
  assign("x_to_evaluate", sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )	, envir = .GlobalEnv )
  
  
  # defining the neighborhood which is to be evaluated	
  # set interaction matrix
  I <- i_matrix	
  # coordinates of neighbours in Interaction matrix I: 
  neighbours_in_I <- which(is.finite(abs(I)/abs(I)), arr.in = TRUE)
  # coordinates relative to the evaluated cell (=  which(is.na(I) ) 
  relrow <- neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
  relcol <- neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
  
  # relative position of the four direct neighbours of a cell
  #interact <- (relrow * dim(X)[2] + relcol)
  
  assign("interact", relrow * dim(X)[2] + relcol, envir = .GlobalEnv )
  
}



####################
## count function ##
####################

count  <- function(x, neighbor) {
  
  neighbors <- numeric(length = prod(x$dim))
  x_logical_with_border <- (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors <- neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}


##################################
## initialise a new "landscape" ##
##################################

init_landscape <- function(states, cover, width = 50, height = 50) {
  if(length(states)!=length(cover)) warning("length of vector 'states' and vector 'cover' differs. I am distributing cover at random!")
  cover <- cover/sum(cover) #normalizing total cover to 1
  cells <- factor(rep(states[1], times = width*height), levels = states) # vector of empty cells in state "0"
  for(i in 2:length(states)) {
    s <- floor(width*height*cover[i]) #how many cells are in state[i]
    cells[sample(which(cells %in% states[1:(i-1)]), s, replace = FALSE)] <- states[i]  # replace s cells, randomly drawn, with state[i]
  }
  
  # wrap landscape object:
  initial <- list(  
    dim = c(width = as.integer(width), height = as.integer(height)),  # first element contains the dimensions of the landscape 
    cells = cells  #contains a random row-wise, factorial vector to fill the grid 
  )
  levels(initial$cells) <- states  #assign cell states 
  class(initial) <- c("list","landscape") # set class of object (required for plotting)
  return(initial)
}


summary.landscape <- function(x) {
  out <- list()
  out$n <- sapply(levels(x$cells), function(y) {sum(x$cells == y)})
  out$cover <- out$n/sum(out$n)
  out$local <- sapply(levels(x$cells), function(y) {mean(count(x,y )/4)})
  return(out)
}

########################################################
## plotting function for objects of class "landscape" ##
########################################################

plot.landscape <- function(x, grid = FALSE, axis = FALSE, cols = "auto", add = FALSE, ani = FALSE, ...) {
  lvls <- levels(x$cells) 
  nlev <- length(lvls)
  if(cols[1] == "auto") cols = c("black", "white")  # default color value
  
  if(ani & Sys.info()[['sysname']] == "Windows") adj = -0.5 else adj = 0 #this adjustment constant is added when producing a pixel accurate png or gif file, requires TRUE when the function is used to plot animated figures. 
  
  if(!add) plot(NA,NA, xlim = c(0.5+adj, x$dim[1]+0.5+adj), ylim = c( x$dim[2]+0.5+adj, 0+0.5+adj), bty = c("n", "o")[grid+1], xaxs = "i", yaxs = "i",xlab = "", ylab = "", xaxt = "n", yaxt = "n", ... ) 
  
  if(axis && !add) axis(3) 
  if(axis && !add) axis(2)
  
  if(grid) border = "grey80" else border = cols[as.numeric(x$cells)]
  
  rect(rep(1:x$dim[1], times = x$dim[2])-.5, rep(1:x$dim[2], each = x$dim[1])-.5, rep(1:x$dim[1], times = x$dim[2])+.5, rep(1:x$dim[2], each = x$dim[1])+.5, col = cols[as.numeric(x$cells)], border = border)
  
  if(grid) box()
}


###################################################
## animated plot for objects of class "caresult" ##
###################################################


animateCA <- function(result, filename) {
  # FIGURE 3 -- animated gif
  library(animation)
  if(Sys.info()[['sysname']] == "Linux") X11.options(antialias = "none") #for Linux Systems to enable pixel-wise plotting in (animated) gif-files. 
  if(Sys.info()[['sysname']] == "Windows") windows.options(antialias = "none") #for Windows Systems to enable pixel-wise plotting in (animated) gif-files. 
  width = result$timeseries[[1]]$dim[1]
  height = result$timeseries[[1]]$dim[2]
  
  saveGIF( 
    for(i in 1:length(result$timeseries) ) {
      par(mar = c(0,0,0,0))
      plot(result$timeseries[[i]], grid = FALSE, ani = TRUE) 
    }
    , movie.name = filename, img.name = "grid", convert = "convert", interval = 0.01/1,
    cmd.fun = system, clean = TRUE, ani.width = width, ani.height = height, outdir = getwd())
}


##################################
## specs for different CA model ##
##################################

# these are updating a landscape object to the next timestep using the model specifications 

### grazing model: Schneider and Kéfi (2015, in review)

update_grazing <- function(x_old, parms_temp, delta = 0.2, subs = 10, timestep = NA) {
  
  x_new <- x_old
  
  for(s in 1:subs) {
    
    
    parms_temp$rho_plus <- sum(x_old$cells == "+")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    parms_temp$Q_plus <- count(x_old, "+")/4  # count local density of occupied fields for each cell:
    
    # 2 - drawing random numbers
    rnum <- runif(x_old$dim[1]*x_old$dim[2]) # one random number between 0 and 1 for each cell
    
    # 3 - setting transition probabilities
    
    recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*1/subs)
    death <- with(parms_temp, (m0+g*(1-protect))*1/subs)
    death[death > 1] <- 1 
    
    regeneration <- with(parms_temp, (r + f*Q_plus)*1/subs)
    
    degradation <- with(parms_temp, (d*1/subs))
    
    # check for sum of probabilities to be inferior 1 and superior 0
    if(any(c(recolonisation+degradation, death, regeneration) > 1 )) warning(paste("a set probability is exceeding 1 in  time step", timestep, "! decrease number of substeps!!!")) 
    if(any(recolonisation < 0)) warning(paste("recolonisation falls below 0 in time step",timestep, "! balance parameters!!!")) 
    if(any(degradation < 0)) warning(paste("degradation falls below 0 in time step",timestep, "! balance parameters!!!"))
    if(any( death < 0)) warning(paste("death falls below 0 in time step",timestep, "! balance parameters!!!")) 
    if(any(regeneration < 0)) warning(paste("regeneration falls below 0 in time step",timestep, "! balance parameters!!!")) 
    
    # 4 - apply transition probabilities  
    
    x_new$cells[which(x_old$cells == "0" & rnum <= recolonisation)] <- "+"
    x_new$cells[which(x_old$cells == "+"  & rnum <= death)] <- "0"
    x_new$cells[which(x_old$cells == "0"  & rnum > recolonisation & rnum <= recolonisation+degradation)] <- "-"  
    x_new$cells[which(x_old$cells == "-"   & rnum <= regeneration)] <- "0"
    
    # 5- store x_new as next x_old
    
    x_old <- x_new
    
  }
  
  ## end of single update call
  return(x_new)
}



################################
## run simulation of CA model ##
################################

ca <- function(x, parms, delta = 0.1, t_max = 1000, t_min = 500, t_eval = 200, isstable = 0.000001, saveeach = 50) {
  
  mapping(x$dim[1], x$dim[2])
  xstats <- summary(x)
  states <- levels(x$cells)
  # calculate timesteps for snapshots:
  n_snaps <- length(seq(t_min-2*t_eval,t_min,saveeach)) # number of snapshots over t_eval
  t_snaps <-  ceiling(t_max/(saveeach*n_snaps))*(n_snaps*saveeach)
  snapshots <- data.frame(time = seq(t_min-2*t_eval, t_max, saveeach), i= seq(t_min-2*t_eval, t_max, saveeach)+1, pos = c(rep(1:n_snaps,2),1) )  
  ### BUG: this is not dynamically adjusting to t_min, t_max, t_eval, saveeach !!!
  
  
  # initialise result object:
  result <- list()  # generate list object
  result$time <- seq(0, t_min) # add vector of realized timesteps
  result$cover <- as.data.frame(t(xstats$cover))
  result$cover <- result$cover[rep(1, t_min+1),] # preallocating memory
  
  result$local <- as.data.frame(t(xstats$local))
  result$local <- result$local[rep(1, t_min+1),] # preallocating memory

  result$timeseries <- list()
  result$timeseries <- lapply(1:n_snaps, function(i) x) # preallocating memory
  
  # --------------- simulation -----------------

  parms_temp <- parms   # set temporary parameter object

  # initialise simulation variables: 
  x_old <- x  # ghost matrix at t_i
  x_new <- x_old
  stability <- 1  # check value for stability
  i = 0  # iterator for simulation timesteps
  
  # starting iterations:
  while(stability > isstable & i <= t_max & result$cover[[1]][i+1] > 0) {

    i <- i +1  # increase iterator
    
    update_grazing(x_old, parms, 10) -> x_new
    
    xstats <- summary(x)
    result$cover[i,] <- xstats$cover
    result$local[i,] <- xstats$local
    
    x_old <- x_new # replace ghost matrix for next iteration

    if(i %in% snapshots$i) {  
      result$timeseries[[match(i, snapshots$i)]] <- x_new
    }
    
    
    if(i > t_min) { # if we are over the minimal timespan 
      t_1 <- (i-2*t_eval):(i-t_eval)-1 # vector of t_eval timesteps previous to the last t_eval timesteps
      t_2 <- (i-t_eval):(i) # vector of the last t_eval timesteps 
      
      if(result$cover[[1]][i] > 0) { 
        stability <- (abs(mean(result$cover[[1]][t_1]) - mean(result$cover[[1]][t_2])))/(mean(result$cover[[1]][t_1])) # calculate stability, i.e. difference between the mean cover in the two evaluation periods t_1 and t_2 
      } else {
        stability <- 0 # set stability to 0 if cover is 0, immediate stop of simulation
      }
      
      result$time[i] <- i # save timestep to results
      
    }
    
  } 
  
  return(result)
}

