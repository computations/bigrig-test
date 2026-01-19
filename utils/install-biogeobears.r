list.of.packages <- c(
  "data.table",
  "ape",
  "rexpokit",
  "cladoRcpp",
  "optimx",
  "phylobase",
  "snow",
  "devtools"
)

.libPaths()

new.packages <- list.of.packages[
  !(list.of.packages %in% installed.packages()[, "Package"])
]
if (length(new.packages)) {
  install.packages(new.packages)
}

library(devtools)
devtools::install_github(
  repo = "nmatzke/BioGeoBEARS",
  build_opts = "--byte-compile"
)
