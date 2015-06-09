README
======

## basic conceptual overview

The code provided here is running spatial disturbance models in a cellular automata framework. The code is written object oriented, meaning that there are objects of class `ca_model` (S3 objects in R) that contain all the model specific information, like the original publication, the update functions, and the cell states. 
These objects can be handled by a function called `ca()` which runs the cellular automata over time. 

Other functions allow the plotting of single snapshots or timeseries, the generation of initial lattices, the calculation of spatial and temporal indicators.

### object structure

A model object is a list `x` that contains

- `x$name`: the name of the model
- `x$ref`: the original reference
- `x$states`: the potential cell states
- `x$cols`: colors for cell states
- `x$parms`: a template for parameters or default parameters 
- `x$update`: the update function


## Bugs

- [x] cut preallocated memory after extinction or set to zero
- [ ] distribution of cells in `init_landscape()`

## ToDo

- [ ] add model selector to `ca()`
  - [ ] implement pred-prey model
  - [x] implement mussel bed model
  - [ ] implement forest gap model
  - [x] check for valid parms object in `ca()`
- [x] automatic color in plot.landscape
- [ ] plot function for results object
- [ ] dynamic checking for equilibrium
- [ ] dynamic calculation of snapshots
- [ ] print.results
- [ ] indicator function
  - [x] select steady state timespan
  - [x] calculate metrics on cover & local cover
  - [ ] calculate metrics on snapshots
- [ ] inter-state clustering