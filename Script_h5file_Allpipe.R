# Batch - Script_h5file_Allpipe.R

#' @title Run all Pipelines

#' @Description This function runs Tapestri and Optima pipelines on same HDF5 file, 
#' converts their outputs to BED format, writes BED files, and then compares the two variant
#' lists for overlaps. The output file in CSV format with common and unique variants for any overlaps
#' 
#'@param HDF5_file
#'@param expected_variants The expected number of variants in the file 
#'@param verbose Logical; if TRUE, prints additional messages during processing.
#'
#'@return A list with elements:
#'\describe{
#'  \item{common_variants}{Charactor vector of variant IDs commonon to both pipelines.}
#'  \item{optima_unique}{charactor vector of variant IDs unique to the Optima pipelin.}
#'  \item{tapestri_unique}{character vector of variant IDs unique to the Tapestri pipeline.}
#'}
#'
#'@Details This function calls internal helper functions to summarise quality metrices, convert
#'variant information to BED format, run each pipeline, and compare the variant lists.

#' @import rhdf5
#' @import dplyr
#' @import VennDiagram
#' @import optima


run_all_pipelines <- function(HDF5_file = NULL, verbose = FALSE){
  
#check if the file exists
if (!file.exists(HDF5_file)) stop("HDF5 file not found.")
  if (verbose) cat("Running Tapestri pipeline...\n")
  tapestri_variants <- run_tapestri_pipeline(HDF5_file, expected_variants)
  write.table(tapestri_variants, file = "passed_1912tapestri_variants.bed",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  if (verbose) cat("Running Optima pipeline...\n")
  optima_variants <- run_optima_pipeline(HDF5_file)
  write.table(optima_variants, file = "passed_1912optima_variants.bed",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  if (verbose) cat("Comparing variant lists...\n")
  comparison <- compare_variant_lists(optima_variants, tapestri_variants)
  
  write.table(data.frame(variant = c(comparison$common_variants, comparison$optima_unique)),
              "optima_all_variants.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(data.frame(variant = c(comparison$common_variants, comparison$tapestri_unique)),
              "tapestri_all_variants.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.csv(data.frame(variant = comparison$common_variants, variant2 = comparison$optima_unique),
            "common_variants.csv", row.names = FALSE)
  write.csv(data.frame(variant = comparison$optima_unique),
            "optima_unique_variants.csv", row.names = FALSE)
  write.csv(data.frame(variant = comparison$tapestri_unique),
            "tapestri_unique_variants.csv", row.names = FALSE)
  
  if (verbose) cat("Pipeline processing and comparison complete.\n")
  return(comparison)
}
  
#' @title Summarize_matrix helper-function
#' 
#'@description A given matrix or vector of quality matrics, the helper function returns a summary
#'vector of length equal to expected number of variants. If the matrix has dimension(cell x variants),
#'it uses colMeans; if transposed(variants x cells), then rowMeans.

#'@param mat A matrix or vector containing quality metrics.
#'@param expected_variants The expected number of variants
  
#'@return A numeric vector of summary values, one per variant.
#'
#'@examples
#'# For a matrix of dimensions 4858 x 82247:
#'# summarised <- summarise_matrix(my_matrix, 82247)
#'
#'@export
  summarize_matrix <- function(mat, expected_variants = 82247) {
  d <- dim(mat)
  if (is.null(d)) {
    # if mat is already a vector then check its length.
    if (length(mat) == expected_variants) {
      return(mat)
    } else {
      stop("Vector length (", length(mat), ") do not match the expected variant count (", expected_variants, ").")
    }
  } else {
    if (d[2] == expected_variants) {
      # Cells x Variants probability, then average over rows.
      return(colMeans(mat, na.rm = TRUE))
    } else if (d[1] == expected_variants) {
      # Variants x Cells probability, then average over columns.
      return(rowMeans(mat, na.rm = TRUE))
    } else {
      stop("Matrix dimensions ", paste(d, collapse = "x"), 
           " do not match expected variant count ", expected_variants)
    }
  }
}

#' BED Format Conversion Function
#' @description Converts data frame having variant information(columns: chr, pos, ref, alt) into
#' BED-format. The start & end coordinates are calculated and constructs an extr_info field.
#' 
#' @param df Data frame having columns \code{Chr}, \code{Pos}, \code{Ref}, and \code{Alt}.
#' 
#' @return Data frame in BED format with columns: \code{Chr}, \code{start_coord}, \code{end_coord}, and \code{extra_info}
#' 
#' @export
convert_to_BED <- function(df) {
  # Here it has been assumed that 'df' has columns as Chr, Pos, Ref, Alt
  bed_df <- df %>%
    mutate(
      start_coord = Pos,
      end_coord = Pos + nchar(Ref) - 1,                         #'end_coord' calculated by adding length of reference allele minus one to 'Pos'
      extra_info = paste("REF=", Ref, "; ALT=", Alt, sep="\t")
    ) %>%
    select(chr = Chr, start_coord, end_coord, extra_info)
  return(bed_df)
}

#' @title Run Tapestri Pipeline
#' Reads HDF5 file, computes per-variant quality metrices and genotypes percentage,
#' filtering threshold is applied, then returns a BED-format data frame of passed variants.
#' @param expected_variants The expected variant numbers
#' 
#' @return BED-format data frame of passed variants
#' 
#' @export
run_tapestri_pipeline <- function(file, expected_variants = 82247) {
  if (!file.exists(file)) stop("Tapestri: h5 file not found.")
  h5closeAll()
  
# Read relevant data from h5 file 
  chrom <- h5read(file, name = "/assays/dna_variants/ca/CHROM")
  pos   <- h5read(file, name = "/assays/dna_variants/ca/POS")
  ref   <- h5read(file, name = "/assays/dna_variants/ca/REF")
  alt   <- h5read(file, name = "/assays/dna_variants/ca/ALT")
  
  # Read relevant data (cells x variants) from h5 file
  DP <- h5read(file, name = "/assays/dna_variants/layers/DP")
  AF <- h5read(file, name = "/assays/dna_variants/layers/AF") * 100
  GQ <- h5read(file, name = "/assays/dna_variants/layers/GQ")
  NGT <- h5read(file, name = "/assays/dna_variants/layers/NGT")
  
  # Replace the genotype value '3' with NA for calculation purposes
  NGT[NGT == 3] <- NA
  
  # Filtering thresholds defined
  DP_cutoff <- 10.           # Minimum read depth required for the variant to considered.
  AF_cutoff <- 20.           # At least 20% of the reads (or alleles) must support the variant
  GQ_cutoff <- 30
  missing_threshold <- 50    # If more than 50% of cells have missing data for that variant, it will be filtered out.
  mutated_threshold <- 10    # At least 10% of cells must be mutated genotype for that variant to pass the filter.

  
# Summarise quality metrices per variant using the helper function
  mean_DP <- summarize_matrix(DP, expected_variants)
  mean_AF <- summarize_matrix(AF, expected_variants)
  mean_GQ <- summarize_matrix(GQ, expected_variants)
  
# Computation of Genotype percentage per variant
# Creates a logical matrix of the same dimension as the NGT matrix
  missing_pct <- summarize_matrix(is.na(NGT), expected_variants) * 100 # Each element is TRUE if corresponding entry in NGT is NA(i.e. missing) & FALSE otherwise
                                                                       # Helper function checks the dimension of its inpu and calculates an average to arrive to the expected number of variants
  mutated_logical <- NGT %in% c(1, 2) # checks each element in the NGT matrix if it is either het(1) or homo(2)
  dim(mutated_logical) <- dim(NGT)  # matrix shape preserved
  mutated_pct <- summarize_matrix(mutated_logical, expected_variants) * 100 # Indicates whether each cell shows a mutated genotype for a given variant & calculates the average (i.e.fraction of cells with a mutated genotype) per variant.
  
# Create a data frame for the variants
  df <- data.frame(
    Chr = chrom,
    Pos = pos,
    Ref = ref,
    Alt = alt,
    DP = mean_DP,
    GQ = mean_GQ,
    AF = mean_AF,
    stringsAsFactors = FALSE
  )
  
  # Identify all variants that meet all the filtering thresholds 
  passed_indices <- which(
    mean_DP >= DP_cutoff &
      mean_AF >= AF_cutoff &
      mean_GQ >= GQ_cutoff &
      missing_pct <= missing_threshold &
      mutated_pct >= mutated_threshold
  )
  
  df <- df[passed_indices, ]
  
  # Convert to BED format
  bed_df <- convert_to_BED(df)
  
  cat("Tapestri pipeline complete. Total passed variants:", nrow(bed_df), "\n")
  return(bed_df)
}

#' @title Run Optima pipeline
#'
#'@discription Read the HDF5 file, filter variants, then converts the filtered variants list into
#'BED format.
#'
#'@return BED-format data frame of filtered variants
#'
#'@export
run_optima_pipeline <- function(file) {
  # Check file exists
  if (!file.exists(file)) stop("Optima: h5 file not found.")
  
  my1912.obj <- readHdf5(
    directory = file,
    sample.name = "1912 Sample",
    omic.type = "DNA"
  )
  
# Filter variants using optima's filterVariant() function
  my1912.obj.filtered <- filterVariant(my1912.obj)
  
# Extract filtered variants 
  filtered_variants <- my1912.obj.filtered@variants
  
# Conversion of filtered variants to a BED-style data frame:computes 'start_coord,end_coord, and extra_info from the existing columns
  variant_df <- data.frame(
    chr = sapply(strsplit(filtered_variants, ":"), `[`, 1),
    start_coord = as.integer(sapply(strsplit(filtered_variants, ":"), `[`, 2)),
    ref_alt = sapply(strsplit(filtered_variants, ":"), `[`, 3),
    stringsAsFactors = FALSE
  )
  
#Split the single column (ref_alt) into separate 'ref' and 'alt' column
  variant_df <- variant_df %>%
    mutate(
      ref = sapply(strsplit(ref_alt, "/"), `[`, 1),
      alt = sapply(strsplit(ref_alt, "/"), `[`, 2),
      end_coord = start_coord + nchar(ref) - 1,
      extra_info = paste("REF=", ref, "; ALT=", alt, sep="\t")
    ) %>%
    select(chr, start_coord, end_coord, extra_info)
  
  cat("Optima pipeline complete. Total passed variants:", nrow(variant_df), "\n")
  return(variant_df)
}

#Variant Lists Function Comparison
compare_variant_lists <- function(optima_df, tapestri_df) {
  
# Creating a unique variant identifiers 
  optima_ids <- with(optima_df, paste(chr, start_coord,  sep=":"))
  tapestri_ids <- with(tapestri_df, paste0("chr",chr, ":",start_coord))
  
#Identify common variants
    common_ids <- intersect(optima_ids, tapestri_ids)

#Identify the unique variants
  optima_unique <- setdiff(optima_ids, tapestri_ids)
  tapestri_unique <- setdiff(tapestri_ids, optima_ids)
  
  cat("Total variants in Optima:", length(optima_ids), "\n")
  cat("Total variants in Tapestri:", length(tapestri_ids), "\n")
  cat("Common variants:", length(common_ids), "\n")
  cat("Unique in Optima:", length(optima_unique), "\n")
  cat("Unique in Tapestri:", length(tapestri_unique), "\n")
  
  # 
  # if (length(optima_ids) > 0 && length(tapestri_ids) > 0) {
  #   venn.diagram(
  #     x = list(Optima = optima_ids, Tapestri = tapestri_ids),
  #     filename = "variant_overlap.svg",
  #     category.names = c("Optima", "Tapestri"),
  #     imagetype = "svg",
  #     height = 800,
  #     width = 800,
  #     resolution = 120,
  #     main = "Venn Diagram of Variant Overlap",
  #     main.cex = 1.5,
  #     main.col = "black",               # Title color
  #     cat.cex = 1.2,                    # Category text size
  #     cex = 1.2,                        # Intersection text size
  #     fill = c("darkorange", "purple"), # Fill colors for the two sets
  #     alpha = 0.5,
  #     cat.col = c("darkblue", "darkgreen")  # Category label colors
  #   )
  # } else {
  #   warning("One or both variant sets are empty; Venn diagram not generated.")
  # }
  # 

#Sends the list as the output of the function 'compare_variant_list().
  return(list(
    common_variants = common_ids,
    optima_unique = optima_unique,
    tapestri_unique = tapestri_unique
  ))
}

# Main Function to run_all_pipelines  
run_all_pipelines <- function() {
  cat("Running Tapestri pipeline...\n") 
  tapestri_variants <- run_tapestri_pipeline(h5_file)
  write.table(tapestri_variants, file = "passed_1912tapestri_variants.bed",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  cat("Running Optima pipeline...\n")  
  optima_variants <- run_optima_pipeline(h5_file)
  write.table(optima_variants, file = "passed_1912optima_variants.bed",
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  cat("Comparing variant lists...\n")
  comparison <- compare_variant_lists(optima_variants, tapestri_variants)
  write.table(data.frame(variant = c(comparison$common_variants, comparison$optima_unique)),
            "optima_all_variants.csv", row.names = FALSE, col.names = FALSE, quote = FALSE) 
  write.table(data.frame(variant = c(comparison$common_variants, comparison$tapestri_unique)),
            "tapestri_all_variants.csv", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.csv(data.frame(variant = comparison$common_variants, variant2 = comparison$optima_unique),
            "common_variants.csv", row.names = FALSE)
  write.csv(data.frame(variant = comparison$optima_unique),
            "optima_unique_variants.csv", row.names = FALSE)
  write.csv(data.frame(variant = comparison$tapestri_unique),
            "tapestri_unique_variants.csv", row.names = FALSE)
  
  cat("Pipeline processing and comparison complete.\n")
  return(comparison)
}
run_all_pipelines()

}
# ----- Run the Entire Pipeline -----
#variant_comparison <- run_all_pipelines()

