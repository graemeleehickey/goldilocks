# Use rmarkdown's vignette styling with a portable, local stylesheet URL.
dir.create("libs", showWarnings = FALSE)
file.copy(
  system.file(
    "rmarkdown/templates/html_vignette/resources/vignette.css",
    package = "rmarkdown",
    mustWork = TRUE
  ),
  "libs/vignette.css",
  overwrite = TRUE
)
