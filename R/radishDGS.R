.onLoad <- function(libname, pkgname) {
  ##this fails with fresh install. Instead, see top of msprime_wrapper()
  #reticulate::source_python(system.file("py/island_model.py",package="radishDGS"),envir=globalenv())
}
#source_python("island_model.py") #NB: mac users may have to modify?

scale_to_0_1 <- function(x) #scales a given vector so that min(x)==0 and max(x)==1
  (x - min(x,na.rm=TRUE))/(max(x,na.rm=TRUE) - min(x,na.rm=TRUE))

lower <- function(x) #returns lower trianglular part of a matrix
  x[lower.tri(x)]

