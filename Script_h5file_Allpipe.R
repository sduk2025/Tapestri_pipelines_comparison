# ===========================================================
# Script: Script_h5file_Allpipe.R
# Description: Runs the Tapestri and Optima pipelines on the same HDF5 file,
#              converts their outputs to BED format, writes the BED files,
#              and then compares the two variant lists for overlaps.
# ===========================================================

# ----- Load Required Libraries -----
library(rhdf5)
library(dplyr)
library(VennDiagram)
library(optima)  # for the Optima pipeline

# ----- Set Working Directory and File Path -----
setwd("/g/data/pq84/single_cell/Tapestri_Sonam/chapter1/Passed_Variants_dir/")
h5_file <- "1912.dna+protein.h5"  # same file for both pipelines
if (!file.exists(h5_file)) stop("HDF5 file not found.")

# ----- Helper Function: summarize_matrix -----
# This function checks the dimensions of a matrix and returns a summary vector
# of length equal to the expected number of variants (82247). If the matrix has
# dimensions [cells, variants], it uses colMeans. If it has dimensions [variants, cells],
# it uses rowMeans.
summarize_matrix <- function(mat, expected_variants = 82247) {
  d <- dim(mat)
  if (is.null(d)) {
    # mat is already a vector; check its length.
    if (length(mat) == expected_variants) {
      return(mat)
    } else {
      stop("Vector length (", length(mat), ") does not match expected variant count (", expected_variants, ").")
    }
  } else {
    if (d[2] == expected_variants) {
      # Likely cells x variants; average over rows.
      return(colMeans(mat, na.rm = TRUE))
    } else if (d[1] == expected_variants) {
      # Likely variants x cells; average over columns.
      return(rowMeans(mat, na.rm = TRUE))
    } else {
      stop("Matrix dimensions ", paste(d, collapse = "x"), 
           " do not match expected variant count ", expected_variants)
    }
  }
}

# ----- Convert to BED Format Function -----
convert_to_BED <- function(df) {
  # Assumes df has columns: Chr, Pos, Ref, Alt
  bed_df <- df %>%
    mutate(
      start_coord = Pos,
      end_coord = Pos + nchar(Ref) - 1,
      extra_info = paste("REF=", Ref, "; ALT=", Alt, sep="\t")
    ) %>%
    select(chr = Chr, start_coord, end_coord, extra_info)
  return(bed_df)
}

# ----- Tapestri Pipeline Function -----
run_tapestri_pipeline <- function(file, expected_variants = 82247) {
  # Check file exists and close any open HDF5 objects
  if (!file.exists(file)) stop("Tapestri: h5 file not found.")
  h5closeAll()
  
  # Read variant metadata from variant-level datasets (should be vectors of length 82247)
  chrom <- h5read(file, name = "/assays/dna_variants/ca/CHROM")
  pos   <- h5read(file, name = "/assays/dna_variants/ca/POS")
  ref   <- h5read(file, name = "/assays/dna_variants/ca/REF")
  alt   <- h5read(file, name = "/assays/dna_variants/ca/ALT")
  
  # Read quality metric matrices (cells x variants)
  DP <- h5read(file, name = "/assays/dna_variants/layers/DP")
  AF <- h5read(file, name = "/assays/dna_variants/layers/AF") * 100
  GQ <- h5read(file, name = "/assays/dna_variants/layers/GQ")
  NGT <- h5read(file, name = "/assays/dna_variants/layers/NGT")
  
  # Replace genotype value '3' with NA for calculation purposes
  NGT[NGT == 3] <- NA
  
  # Define filtering thresholds
  DP_cutoff <- 10
  AF_cutoff <- 20
  GQ_cutoff <- 30
  missing_threshold <- 50    # maximum allowed % missing
  mutated_threshold <- 10    # minimum required % mutated
  
  # Summarize quality metrics per variant using the helper function:
  mean_DP <- summarize_matrix(DP, expected_variants)
  mean_AF <- summarize_matrix(AF, expected_variants)
  mean_GQ <- summarize_matrix(GQ, expected_variants)
  
  # Genotype percentages per variant Computation:
  missing_pct <- summarize_matrix(is.na(NGT), expected_variants) * 100
  
  mutated_logical <- NGT %in% c(1, 2)
  dim(mutated_logical) <- dim(NGT)  # Preserve the matrix shape
  mutated_pct <- summarize_matrix(mutated_logical, expected_variants) * 100
  
  # A data frame with variant-level data build
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
  
  # Applying quality and genotype filters
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

# ----- Optima Pipeline Function -----#
run_optima_pipeline <- function(file) {
  # Check file exists
  if (!file.exists(file)) stop("Optima: h5 file not found.")
  
  # Optima package to read in the h5 file
  my1912.obj <- readHdf5(
    directory = file,
    sample.name = "1912 Sample",
    omic.type = "DNA"
  )
  
  # Filter variants using optima's filterVariant() function
  my1912.obj.filtered <- filterVariant(my1912.obj)
  
  # Extract filtered variants (assumed format: "chr:pos:ref/alt")
  filtered_variants <- my1912.obj.filtered@variants
  
  # Convert filtered variants to a BED-style data frame:
  variant_df <- data.frame(
    chr = sapply(strsplit(filtered_variants, ":"), `[`, 1),
    start_coord = as.integer(sapply(strsplit(filtered_variants, ":"), `[`, 2)),
    ref_alt = sapply(strsplit(filtered_variants, ":"), `[`, 3),
    stringsAsFactors = FALSE
  )
  
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

# ----- Variant Lists Function Comparison-----
compare_variant_lists <- function(optima_df, tapestri_df) {
  # Create unique variant identifiers (e.g., "chr:start:extra_info")
  optima_ids <- with(optima_df, paste(chr, start_coord,  sep=":"))
  tapestri_ids <- with(tapestri_df, paste0("chr",chr, ":",start_coord))
  
  common_ids <- intersect(optima_ids, tapestri_ids)
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
  return(list(
    common_variants = common_ids,
    optima_unique = optima_unique,
    tapestri_unique = tapestri_unique
  ))
}

# ----- Main Function: run_all_pipelines -----
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

# ----- Run the Entire Pipeline -----
variant_comparison <- run_all_pipelines()

