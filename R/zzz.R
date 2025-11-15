
.onLoad <- function(libname, pkgname) {
  rp <- Sys.getenv("RETICULATE_PYTHON")
  if (nzchar(rp)) {
    reticulate::use_python(rp, required = TRUE)
  } else {
    print("RETICULATE_PYTHON not set, using default Python environment")
  }
}
