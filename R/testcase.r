source("R/functions.r")
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

