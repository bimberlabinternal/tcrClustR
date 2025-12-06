
.onLoad <- function(libname, pkgname) {
  rp <- Sys.getenv("RETICULATE_PYTHON")
  if (nzchar(rp)) {
    reticulate::use_python(rp, required = TRUE)
  } else {
    packageStartupMessage("RETICULATE_PYTHON not set, using default Python environment")
  }
}

#' Setup Python Environment for tcrClustR
#'
#' @description
#' This function validates the Python installation, checks for required modules,
#' and installs dependencies if needed. It's designed to be called before any
#' tcrClustR functions that rely on Python (e.g., RunTcrdist3, .FormatMetadata).
#'
#' @param pythonExecutable Path to the Python executable. If NULL (default), will
#'   use reticulate::py_exe() or fall back to system Python3.
#' @param installMissing Logical. If TRUE (default), will attempt to install missing
#'   Python packages. If FALSE, will only validate and report missing packages.
#' @param useRequirementsTxt Logical. If TRUE (default), will attempt to install
#'   dependencies from requirements.txt if available. If FALSE, will install
#'   individual packages.
#' @param verbose Logical. If TRUE, will display detailed progress messages.
#'
#' @return Invisible logical: TRUE if all dependencies are available, FALSE otherwise.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Validates Python executable exists and is functional
#'   \item Checks for required Python modules: tcrdist3, pandas, numpy, rpy2
#'   \item If installMissing=TRUE, installs missing packages via pip
#'   \item Returns status of dependency availability
#' }
#'
#' The function uses requirements.txt for dependency management when available,
#' and falls back to individual package installation if not found.
#'
#' @examples
#' \dontrun{
#' # Check and install Python dependencies
#' SetupPythonEnvironment()
#'
#' # Just validate without installing
#' SetupPythonEnvironment(installMissing = FALSE)
#'
#' # Use specific Python executable
#' SetupPythonEnvironment(pythonExecutable = "/usr/bin/python3")
#' }
#'
#' @export
SetupPythonEnvironment <- function(pythonExecutable = NULL,
                                    installMissing = TRUE,
                                    useRequirementsTxt = TRUE,
                                    verbose = TRUE) {
  #determine python executable
  if (is.null(pythonExecutable)) {
    pythonExecutable <- tryCatch({
      reticulate::py_exe()
    }, error = function(e) NULL)
    
    if (is.null(pythonExecutable) || pythonExecutable == "") {
      pythonExecutable <- Sys.which("python3")
      if (pythonExecutable == "") {
        pythonExecutable <- Sys.which("python")
      }
    }
  }
  
  if (pythonExecutable == "" || !file.exists(pythonExecutable)) {
    stop("No valid Python executable found. Please install Python 3.8+ or specify pythonExecutable parameter.")
  }
  
  if (verbose) {
    message("Using Python executable: ", pythonExecutable)
  }
  
  #check python version
  version_output <- tryCatch({
    system2(pythonExecutable, "--version", stdout = TRUE, stderr = TRUE)
  }, error = function(e) {
    stop("Failed to execute Python. Error: ", e$message)
  })
  
  if (verbose) {
    message("Python version: ", paste(version_output, collapse = " "))
  }
  
  # NOTE: The package is installed as 'tcrdist3' but imported as 'tcrdist'
  required_modules <- c("tcrdist", "pandas", "numpy", "rpy2")
  
  #check which modules are missing
  missing_modules <- character(0)
  for (mod in required_modules) {
    result <- system2(pythonExecutable, 
                     c("-c", shQuote(paste0("import ", mod))),
                     stdout = FALSE, 
                     stderr = FALSE)
    
    if (result != 0) {
      missing_modules <- c(missing_modules, mod)
    }
  }
  
  if (length(missing_modules) > 0) {
    if (verbose) {
      message("Missing Python modules: ", paste(missing_modules, collapse = ", "))
    }
    
    if (installMissing) {
      if (verbose) {
        message("Installing missing Python dependencies...")
      }
      
      if (useRequirementsTxt) {
        #locate requirements.txt
        requirements_path <- system.file("requirements.txt", package = "tcrClustR")
        
        #if not found in installed package, look in source directory
        if (requirements_path == "" || !file.exists(requirements_path)) {
          pkg_root <- system.file(package = "tcrClustR")
          if (pkg_root != "") {
            requirements_path <- file.path(dirname(pkg_root), "requirements.txt")
          }
          
          if (!file.exists(requirements_path)) {
            inst_dir <- system.file("inst", package = "tcrClustR")
            if (inst_dir != "") {
              requirements_path <- file.path(dirname(inst_dir), "requirements.txt")
            }
          }
        }
        
        if (file.exists(requirements_path)) {
          if (verbose) {
            message("Installing from requirements.txt: ", requirements_path)
          }
          
          #upgrade pip
          install_result <- system2(
            pythonExecutable,
            c("-m", "pip", "install", "--upgrade", "pip"),
            stdout = TRUE,
            stderr = TRUE
          )
          
          if (verbose) {
            cat(install_result, sep = "\n")
          }
          
          # Install from requirements.txt
          install_result <- system2(
            pythonExecutable,
            c("-m", "pip", "install", "-r", shQuote(requirements_path)),
            stdout = TRUE,
            stderr = TRUE
          )
          
          if (verbose) {
            cat(install_result, sep = "\n")
          }
        } else {
          warning("requirements.txt not found. Falling back to individual package installation.")
          useRequirementsTxt <- FALSE
        }
      }
      
      if (!useRequirementsTxt) {
        #install individual packages
        packages_to_install <- c()
        if ("tcrdist3" %in% missing_modules) {
          packages_to_install <- c(packages_to_install, "tcrdist3>=0.2.2")
        }
        if ("pandas" %in% missing_modules) {
          packages_to_install <- c(packages_to_install, "pandas>=1.1.0")
        }
        if ("numpy" %in% missing_modules) {
          packages_to_install <- c(packages_to_install, "numpy>=1.19.0")
        }
        if ("rpy2" %in% missing_modules) {
          packages_to_install <- c(packages_to_install, "rpy2>=3.4.0")
        }
        
        for (pkg in packages_to_install) {
          if (verbose) {
            message("Installing ", pkg, "...")
          }
          
          install_result <- system2(
            pythonExecutable,
            c("-m", "pip", "install", shQuote(pkg)),
            stdout = TRUE,
            stderr = TRUE
          )
          
          if (verbose) {
            cat(install_result, sep = "\n")
          }
        }
      }
      
      #re-check modules after installation
      still_missing <- character(0)
      for (mod in required_modules) {
        result <- system2(pythonExecutable, 
                         c("-c", shQuote(paste0("import ", mod))),
                         stdout = FALSE, 
                         stderr = FALSE)
        if (result != 0) {
          still_missing <- c(still_missing, mod)
        }
      }
      
      if (length(still_missing) > 0) {
        warning("Failed to install some Python modules: ", paste(still_missing, collapse = ", "))
        return(invisible(FALSE))
      }
      
      if (verbose) {
        message("All Python dependencies successfully installed!")
      }
    } else {
      warning("Missing required Python modules: ", paste(missing_modules, collapse = ", "),
              "\nRun SetupPythonEnvironment(installMissing = TRUE) to install them.")
      return(invisible(FALSE))
    }
  } else {
    if (verbose) {
      message("All required Python modules are available.")
    }
  }
  
  return(invisible(TRUE))
}
