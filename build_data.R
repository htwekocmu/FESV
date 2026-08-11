################################################################################
# build_data.R
#
# Builds the example datasets shipped with the FESV package from the raw CSVs
# of the empirical application (World Bank WDI, 2000-2023).
#
# Run from the package root with:
#   Rscript data-raw/build_data.R
################################################################################

make_panel <- function(path, tag) {
  dat <- read.csv(path)
  dat <- dat[order(dat$id, dat$year), ]
  keep <- c("id", "year", "Income", "FDI", "INF", "POP", "TOP", "UR", "EXR")
  dat <- dat[, intersect(keep, colnames(dat))]
  cat(tag, ": N =", length(unique(dat$id)),
      "| T =", length(unique(dat$year)),
      "| NT =", nrow(dat), "\n")
  dat
}

developing_panel <- make_panel("data-raw/developing.csv", "developing_panel")
developed_panel  <- make_panel("data-raw/developed.csv",  "developed_panel")
all_panel        <- make_panel("data-raw/all.csv",        "all_panel")

save(developing_panel, file = "data/developing_panel.rda", compress = "bzip2")
save(developed_panel,  file = "data/developed_panel.rda",  compress = "bzip2")
save(all_panel,        file = "data/all_panel.rda",        compress = "bzip2")
cat("saved data/*.rda\n")
