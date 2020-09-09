.onLoad <- function(libname, pkgname) {
  ##this fails with fresh install. Instead, see top of msprime_wrapper()
  #reticulate::source_python(system.file("py/island_model.py",package="radishDGS"),envir=globalenv())
}
#source_python("island_model.py") #NB: mac users may have to modify?
#' @export


