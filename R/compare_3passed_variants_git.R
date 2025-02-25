#Batch - Compare Filtered Variant Files
# 
#'@discription This function takes in input variants from the Tapestri, Optima and Inhouse pipelines,
#' compares the three variant lists for any overlaps. It generates Venndiagram and Upset plot and writes
#' out text files - common, unique to each of the pipeline.
#' 
#' @param tapestri_file Character string with filename for Tapestri variant list
#' \description{
#'    \item{common_variants} {A character vector of variant keys common to all pipeline.
#'    \item{tapestri_unique}{A character vector of variant keys unique to the Tapestri pipeline.}
#'    \item{optima_unique}{A character vector of variant keys unique to the Optima pipeline.}
#'    \item{inhouse_unique}{A character vector of variant keys unique to the Inhouse pipeline}
#'}
#'
#' @examples {Rscript compare_3passed_variants.R <Tapestri_keys.txt> <Optima_keys.txt>\n"}
#' 
#' @export

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(VennDiagram)   
  library(UpSetR)        
})

compare_3passed_variants <- function(optima_variants = NULL, tapestri_variants = NULL, inhouse_variants, verbose = FALSE){

  # Check if files exist or not
  if (!file.exists(fileTapestri)) stop(paste("File not found:", fileTapestri))
  if (!file.exists(fileOptima)) stop(paste("File not found:", fileOptima))
  if (!file.exists(fileInhouse)) stop(paste("File not found:", fileInhouse))
  
  cat("Reading variant key files...\n")
  setTapestri <- unique(read_lines(fileTapestri))
  setOptima   <- unique(read_lines(fileOptima))
  setInhouse <- unique(read_lines(fileInhouse))
  
  cat("File sizes:\n")
  cat("  Tapestri:", length(setTapestri), "variants\n")
  cat("  Optima:", length(setOptima), "variants\n")
  cat(" Inhouse:", length(setInhouse))
  
  #Create named list for VennDiagram & UpSet plot 
  variantList <- list(
    Tapestri = setTapestri,
    Optima   = setOptima,
    Inhouse = setInhouse
  )
  
cat("File sizes:\n")
cat("  Tapestri:", length(setTapestri), "variants\n")
cat("  Optima:", length(setOptima), "variants\n")
cat(" Inhouse:", length(setInhouse))

#Create named list for VennDiagram & UpSet 
variantList <- list(
  Tapestri = setTapestri,
  Optima   = setOptima,
  Inhouse = setInhouse
)

# Venn Diagram Generation
cat("Generating Venn diagram...\n")
vennFile <- "v_diagram_3pipe_variants.png"
venn.diagram(
  x = variantList,
  filename = vennFile,
  imagetype = "png",
  height = 800,
  width  = 800,
  resolution = 120,
  main = "Venn Diagram: 3pipelines",
  main.cex = 1.5,
  cat.cex = 1.2,
  cex = 1.2,
  fill = c("blue", "darkblue", "darkgreen"),  # Three colors for three sets
  alpha = 0.5,
  cat.col = c("cyan", "black", "red")
)
cat("Venn diagram saved to:", vennFile, "\n")

# Generate UpSet plot
cat("Generating UpSet plot...\n")
upsetFile <- "upset_3pipes_variants.pdf"
pdf(upsetFile, width = 7, height = 5)
mat <- fromList(variantList)
upset(
  mat,
  sets = c("Tapestri", "Optima", "Inhouse"),
  order.by = "freq",
  main.bar.color = "darkblue",
  sets.bar.color = "lightblue",
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
# Intersection of all 3 pipes:
intersect_3 <- intersect(setTapestri, intersect(setOptima, setInhouse))
intersect_Tap_Opt <- setdiff(intersect(setTapestri, setOptima), setInhouse)
intersect_Tap_Inh <- setdiff(intersect(setTapestri, setInhouse), setOptima)
intersect_Opt_Inh <- setdiff(intersect(setOptima, setInhouse), setTapestri)
unique_Tap <- setdiff(setTapestri, union(setOptima, setInhouse))
unique_Opt <- setdiff(setOptima, union(setTapestri, setInhouse))
unique_Inh <- setdiff(setInhouse, union(setTapestri, setOptima))

cat("Writing subset lists to files...\n")
write_lines(intersect_3,       "venn_subset_all3.txt")      
write_lines(intersect_Tap_Opt, "venn_subset_TapestriOptima.txt")  
write_lines(intersect_Opt_Inh, "venn_subset_OptimaInhouse.txt")   
write_lines(intersect_Tap_Inh, "venn_subset_TapestriInhouse.txt") 
write_lines(unique_Opt,        "venn_subset_OptimaOnly.txt")      
write_lines(unique_Inh,        "venn_subset_InhouseOnly.txt")     
write_lines(unique_Tap,        "venn_subset_TapestriOnly.txt")    

cat("\nDone! You now have:\n",
    "  1) venn_3algos.png\n",
    "  2) upset_3algos.pdf\n",
    "  3) The 7 subset files:\n",
    "     venn_subset_all3.txt\n",
    "     venn_subset_TapestriOptima.txt\n",
    "     venn_subset_OptimaInhouse.txt\n",
    "     venn_subset_TapestriInhouse.txt\n",
    "     venn_subset_OptimaOnly.txt\n",
    "     venn_subset_InhouseOnly.txt\n",
    "     venn_subset_TapestriOnly.txt\n",
    sep="")

return(list(
  intersect_3 = common_variants,
  tapestri_unique = unique_tapestri,
  optima_unique = unique_Optima,
  inhouse_unique = unique_Inhouse
  ))
}

#Run the script only if it is executed from the command lin
#Check whether the number of command-line arguments is exactly 3
if (length(args) != 3) {
  args <- commandArgs(trailingOnly = TRUE) # returns a vector which is stored in 'args'.
  if (length(args) != 3) {
  cat("Usage: Rscript compare_3passed_variants_git.R <Tapestri_keys.txt> <Optima_keys.txt> <Inhouse_keys.txt>\n")
  quit(status = 1) #script exiting with an error code stopping the script from running further
  }
}

#Take command-line arguments provided by the user.
fileTapestri <- args[1]  # e.g. "Tapestri_passed1912_variants.txt"
fileOptima   <- args[2]  # e.g. "Optima_passed1912_variants.txt"
fileInhouse <- args[3]   # e.g. "Inhouse_passed1912_variants.txt"

comparison <- compare_variant_files(fileTapestri,fileOptima, fileInhouse)

  