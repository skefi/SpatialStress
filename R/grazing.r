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


### grazing model: Schneider and Kéfi (2015, in review)

grazing <- list()
class(grazing) <- "ca_model"
grazing$name <- "Spatial Grazing Model"
grazing$ref <- "Schneider and Kéfi 2015, in review"
grazing$states <- c("+", "0", "-")
grazing$cols <- grayscale(3)
grazing$parms <- list(
  del = 0.9, # local seed dispersal
  b = 0.2, # environmental quality
  c_ = 0.2, # global competition
  m0 = 0.05, # intrinsic mortality
  g = 0.2, # grazing pressure
  r = 0.01, # regeneration rate of degraded cells
  f = 0.9, # local facilitation
  d = 0.2, # intrinsic degradation rate 
  protect = 0.9 # associational resistance against grazing
) 
grazing$update <- function(x_old, parms_temp, delta = 0.2, subs = 10, timestep = NA) {
  
  x_new <- x_old
  
  for(s in 1:subs) {
    
    
    parms_temp$rho_plus <- sum(x_old$cells == "+")/(x_old$dim[1]*x_old$dim[2]) # get initial vegetation cover
    parms_temp$Q_plus <- count(x_old, "+")/4  # count local density of occupied fields for each cell:
    
    # 2 - drawing random numbers
    rnum <- runif(x_old$dim[1]*x_old$dim[2]) # one random number between 0 and 1 for each cell
    
    # 3 - setting transition probabilities
    
    if(parms_temp$rho_plus > 0) {
      recolonisation <- with(parms_temp, (del*rho_plus+(1-del)*Q_plus)*(b-c_*rho_plus)*1/subs)
      
      death <- with(parms_temp, (m0+g*(1-protect))*1/subs)
      death[death > 1] <- 1 
    } else {
      recolonisation <- 0
      death <- 1
    }
    
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


