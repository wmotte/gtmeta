# Build the standalone IGMI-styled HTML site from the package vignette.
# Run from the package root:  Rscript docs/build.R
root <- normalizePath(".")
stopifnot(file.exists(file.path(root, "DESCRIPTION")))

rmarkdown::render(
  file.path(root, "vignettes", "igmi-tour.Rmd"),
  output_format = rmarkdown::html_document(
    theme = NULL,
    highlight = "pygments",
    css = file.path(root, "docs", "igmi.css"),
    self_contained = TRUE,
    toc = TRUE,
    toc_depth = 2,
    includes = rmarkdown::includes(
      before_body = file.path(root, "docs", "header.html"),
      after_body  = file.path(root, "docs", "footer.html")
    )
  ),
  output_file = "index.html",
  output_dir = file.path(root, "docs"),
  knit_root_dir = root
)
