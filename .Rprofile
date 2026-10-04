# Keep the project library and renv state local to this checkout. This also
# avoids inheriting a user-wide library location during project activation.
Sys.setenv(
  RENV_PATHS_LIBRARY = file.path(getwd(), "renv", "library"),
  RENV_PATHS_ROOT = file.path(getwd(), "renv", ".state"),
  # Use R's standard library directly rather than copying it into a sandbox.
  # Project dependencies still resolve from the private project library first.
  RENV_CONFIG_SANDBOX_ENABLED = "FALSE"
)
source("renv/activate.R")
