#NOTE: This file runs before any tests are executed in testthat.
#It configures the Python environment for all test files.
#This is primarily for CI/CD compatibility (GitHub Actions) but also works locally.
#Users running tests locally do not need to manually source this file.

#get RETICULATE_PYTHON from environment (set by GitHub Actions workflow or local .Renviron)
reticulate_python <- Sys.getenv("RETICULATE_PYTHON")

#if not set, try to auto-detect python3 or python
if (!nzchar(reticulate_python)) {
  reticulate_python <- Sys.which("python3")
  if (!nzchar(reticulate_python)) {
    reticulate_python <- Sys.which("python")
  }
}

#set RETICULATE_PYTHON for all tests
if (nzchar(reticulate_python)) {
  Sys.setenv(RETICULATE_PYTHON = reticulate_python)
  message("Test setup: RETICULATE_PYTHON set to ", reticulate_python)
  
  #verify Python is accessible
  tryCatch({
    version <- system2(reticulate_python, "--version", stdout = TRUE, stderr = TRUE)
    message("Test setup: Python version: ", paste(version, collapse = " "))
  }, error = function(e) {
    warning("Test setup: Failed to verify Python: ", e$message)
  })
} else {
  message("Test setup: No Python executable found - Python-dependent tests will be skipped")
}
