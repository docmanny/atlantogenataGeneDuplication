if (!requireNamespace("checkpoint", quietly = TRUE)) {install.packages("checkpoint")}
if (!requireNamespace("devtools", quietly = TRUE)) {install.packages("devtools")}

library(checkpoint)
library(devtools)

dir.create(file.path(tempdir(), ".checkpoint"), recursive = TRUE, showWarnings = FALSE)

# Note: PythonInR fails if compiling from CRAN, use bitbucket version instead
Sys.setenv(TAR="/bin/tar")
install_bitbucket("Floooo/PythonInR@3d88749a2d22ec4f7479c9802ffccb86f81f9538")

# Get checkpoint
checkpoint("2019-01-04", checkpointLocation = tempdir(), )

CRAN.packages <- c(
  "tidyverse",
  "ape",
  "reshape2",
  "workflowr",
  "WebGestaltR",
  "UpSetR",
  "scales",
  "ggpubr",
  "ggthemes",
  "ggimage",
  "RColorBrewer",
  "digest",
  "reticulate",
  "BiocManager"
  )

# Make sure all packages were really installed by checkpoint()
for (pkg in CRAN.packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {install.packages(pkg)}
}

bc.packages <- c(
  "ggtree",
  "GenomicRanges",
  "GenomicFeatures"
)

BiocManager::install(ask=F, update=F, version="3.8")

#BiocInstaller::biocinstallRepos()


