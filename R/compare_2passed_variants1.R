
# Batch - Compare Passed Variant Files

# 
#' @description This function takes input variants list from the Tapestri and Optima pipelines,
#' compares the two variant lists for any overlaps.It generates Venndiagram and Upset Plot and writes
#' out text files - common, unique to each pipeline.
#'
#'@param tapestri_file Character string with the filename for Tapestri variant list 
#'e.g. "1912_Optima_keys.txt and 1912_Tapestri_keys.txt"

#'@return List with elements:
#'\describe{
#'  \item{common_variants}{A character vector of variant keys common to both pipelines.}
#'  \item{tapestri_unique}{A character vector of variant keys unique to the Tapestri pipeline.}
#'  \item{optima_unique}{A character vector of variant keys unique to the Optima pipeline.}
#'}
#'
#' @examples {Rscript compare_passed_variants.R <Tapestri_keys.txt> <Optima_keys.txt>\n"}
#' 
#' @export

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(VennDiagram)   
  library(UpSetR)        
})

compare_2passed_variants <- function(optima_variants = NULL, tapestri_variants = NULL, verbose = FALSE){

# Check if files exist or not
  if (!file.exists(fileTapestri)) stop(paste("File not found:", fileTapestri))
  if (!file.exists(fileOptima))   stop(paste("File not found:", fileOptima))

  cat("Reading variant key files...\n")
  setTapestri <- unique(read_lines(fileTapestri))
  setOptima   <- unique(read_lines(fileOptima))
  
  cat("File sizes:\n")
  cat("  Tapestri:", length(setTapestri), "variants\n")
  cat("  Optima:   ", length(setOptima), "variants\n")
  
#Create named list for VennDiagram & UpSet 
  variantList <- list(
    Tapestri = setTapestri,
    Optima   = setOptima
  )
  
# Venn Diagram Generation
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

# Identify variant subsets
common_variants <- intersect(setTapestri, setOptima)
unique_Tapestri <- setdiff(setTapestri, setOptima)
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

return(list(
  common_variants = common_variants,
  tapestri_unique = unique_Tapestri,
  optima_unique = unique_Optima
))
}

# Run the script only if it is executed from the command line
#Checks whether the number of command-line arguments is exactly 2
if (length(args) != 2) {
  args <- commandArgs(trailingOnly = TRUE) 
  if (length(args) != 2) {
    cat("Usage: Rscript compare_2passed_variants.R <Tapestri_keys.txt> <Optima_keys.txt>\n")
    quit(status = 1)
  }
  cat("Usage: Rscript compare_passed_variants.R <Tapestri_keys.txt> <Optima_keys.txt>\n")
  quit(status = 1) #script exiting with an error code stopping the script from running further
}
#Take command-line arguments provided by the user.
fileTapestri <- args[1]  # e.g. "1912_Tapestri_keys.txt"
fileOptima   <- args[2]  # e.g. "1912_Optima_keys.txt"

comparison <- compare_variant_files(fileTapestri,fileOptima)
