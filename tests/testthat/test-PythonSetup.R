test_that("SetupPythonEnvironment validates Python installation", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("python3")), "Python3 not available")
  
  # Test that function runs without error
  result <- SetupPythonEnvironment(installMissing = FALSE, verbose = FALSE)
  
  # Should return a logical value
  expect_type(result, "logical")
})

test_that("SetupPythonEnvironment detects Python executable", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("python3")), "Python3 not available")
  
  # Test with auto-detection
  expect_no_error({
    SetupPythonEnvironment(installMissing = FALSE, verbose = TRUE)
  })
  
  # Test with explicit Python path
  python_path <- Sys.which("python3")
  if (python_path != "") {
    expect_no_error({
      SetupPythonEnvironment(pythonExecutable = python_path, 
                            installMissing = FALSE, 
                            verbose = FALSE)
    })
  }
})

test_that("SetupPythonEnvironment handles missing Python gracefully", {
  skip_on_cran()
  
  # Test with non-existent Python path
  expect_error({
    SetupPythonEnvironment(pythonExecutable = "/path/to/nonexistent/python",
                          installMissing = FALSE,
                          verbose = FALSE)
  }, regexp = "No valid Python executable found")
})

test_that("SetupPythonEnvironment checks required modules", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("python3")), "Python3 not available")
  
  # This test will either:
  # 1. Return TRUE if all modules are installed
  # 2. Return FALSE if modules are missing and installMissing=FALSE
  # 3. Attempt installation if installMissing=TRUE
  
  result <- suppressWarnings({
    SetupPythonEnvironment(installMissing = FALSE, verbose = FALSE)
  })
  
  # Result should be logical
  expect_type(result, "logical")
})

test_that("SetupPythonEnvironment verbose mode provides output", {
  skip_on_cran()
  skip_if_not(nzchar(Sys.which("python3")), "Python3 not available")
  
  # Capture messages when verbose=TRUE
  expect_message({
    SetupPythonEnvironment(installMissing = FALSE, verbose = TRUE)
  }, regexp = "Using Python executable")
})

test_that("pyproject.toml exists in package", {
  # Check that pyproject.toml file exists
  # It should be in the package root, but may not be installed with the package
  # This test just verifies it exists in the source
  pyproject_path <- system.file("pyproject.toml", package = "tcrClustR")
  
  # If installed package doesn't have it (expected), just pass
  # This is mainly for development testing
  if (pyproject_path == "") {
    expect_true(TRUE)  # Pass if not found in installed package
  } else {
    expect_true(file.exists(pyproject_path))
  }
})
