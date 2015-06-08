
inittest <- init_landscape(c("+","0","-"), c(0.5,0.4,0.1))
mapping(50,50)
summary(inittest)

parmstest <- list(
del = 0.9,
b = 0.8,
c_ = 0.2,
m0 = 0.05,
g = 0.2,
r = 0.01,
f = 0.9,
d = 0.1,
protect = 0.5

  ) 

simtest <- ca(inittest, parmstest)


update_grazing <- function(x_old, parms_temp, delta = 0.2, subs = 10, timestep = NA) {
  
  x_new <- x_old
  
  for(s in 1:subs) {
    
    
    parms_temp$rho_plus <- sum(x_old$cells == "+")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    parms_temp$Q_plus <- count(x_old, "+")/4  # count local density of occupied fields for each cell:
    
    # 2 - drawing random numbers
    rnum <- runif(width*height) # one random number between 0 and 1 for each cell
    
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
  
update_grazing(inittest, parmstest, 10) -> x_new
  