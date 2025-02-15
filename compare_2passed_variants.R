
# ===========================================================
# Script: compare_2passed_variants.R
# Description: Takes in input variants from the Tapestri and Optima pipelines,
#              compares the two variant lists for any overlaps.
# ===========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(VennDiagram)   # for venn.diagram
  library(UpSetR)        # for upset() and fromList()
})

# Set working directory (adjust as needed)
setwd("/g/data/pq84/single_cell/Tapestri_Sonam/chapter1/Passed_Variants_dir/")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript compare_passed_variants.R <Tapestri_keys.txt> <Optima_keys.txt>\n")
  quit(status = 1)
}

fileTapestri <- args[1]  # e.g. "1912_Tapestri_keys.txt"
fileOptima   <- args[2]  # e.g. "1912_Optima_keys.txt"

# Check if files exist
if (!file.exists(fileTapestri)) stop(paste("File not found:", fileTapestri))
if (!file.exists(fileOptima))   stop(paste("File not found:", fileOptima))

cat("Reading variant key files...\n")
setTapestri <- unique(read_lines(fileTapestri))
setOptima   <- unique(read_lines(fileOptima))

cat("File sizes:\n")
cat("  Tapestri:", length(setTapestri), "variants\n")
cat("  Optima:   ", length(setOptima), "variants\n")

# Create a named list for VennDiagram & UpSet
variantList <- list(
  Tapestri = setTapestri,
  Optima   = setOptima
)

# Generate Venn Diagram
cat("Generating Venn diagram...\n")
vennFile <- "venn_2algos_passed.png"
venn.diagram(
  x = variantList,
  filename = vennFile,
  imagetype = "png",
  height = 800,
  width  = 800,
  resolution = 120,
  main = "Venn Diagram of Passed Variant Lists",
  main.cex = 1.5,
  cat.cex = 1.2,
  cex = 1.2,
  fill = c("purple", "#377EB8"),  # Two colors for two sets
  alpha = 0.5,
  cat.col = c("darkgreen", "#377EB8")
)
cat("Venn diagram saved to:", vennFile, "\n")

# Generate UpSet plot
cat("Generating UpSet plot...\n")
upsetFile <- "upset_2algos.pdf"
pdf(upsetFile, width = 7, height = 5)
mat <- fromList(variantList)
upset(
  mat,
  sets = c("Tapestri", "Optima"),
  order.by = "freq",
  main.bar.color = "darkred",
  sets.bar.color = "steelblue",
  number.angles = 0,
  text.scale = 1.4,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  mb.ratio = c(0.6, 0.4),
  point.size = 3.5,
  line.size = 1.0
)
dev.off()
cat("UpSet plot saved to:", upsetFile, "\n")

# Identify each Venn subset
# Common to both pipelines
common_variants <- intersect(setTapestri, setOptima)

# Unique to Tapestri only
unique_Tapestri <- setdiff(setTapestri, setOptima)

# Unique to Optima only
unique_Optima <- setdiff(setOptima, setTapestri)

cat("Writing subset lists to files...\n")
write_lines(common_variants, "venn_subset_common.txt")
write_lines(unique_Tapestri, "venn_subset_TapestriOnly.txt")
write_lines(unique_Optima, "venn_subset_OptimaOnly.txt")

cat("\nDone! The following files have been created:\n",
    "  1) ", vennFile, "\n",
    "  2) ", upsetFile, "\n",
    "  3) venn_subset_common.txt\n",
    "  4) venn_subset_TapestriOnly.txt\n",
    "  5) venn_subset_OptimaOnly.txt\n",
    sep = "")
