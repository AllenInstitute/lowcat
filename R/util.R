#' Installs bioconductor dependencies
install_bioc_deps <- function() {
  if(!"BiocManager" %in% installed.packages()) {
    install.packages("BiocManager")
  }

  bioc_packages <- c("GenomicAlignments",
                     "GenomicRanges",
                     "rsamtools",
                     "rtracklayer")

  missing_bioc_packages <- setdiff(bioc_packages, installed.packages())

  BiocManager::install(missing_bioc_packages)
}
