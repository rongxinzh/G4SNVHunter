#' Load Small Variant Data from VCF or MAF Files
#'
#' This function loads variant data from either a standard VCF or MAF file and
#' filters for small variants, such as SNVs, INDELs, and DELINs.
#'
#' @param variant_file A \code{character} string specifying the path to the
#' variant file, which can be in either VCF or MAF format.
#' @param file_type A \code{character} string specifying the type of the
#' input file: either \code{"vcf"} or \code{"maf"}. This parameter must be
#' provided. If \code{"maf"} is specified, the following columns must be
#' present in the MAF file: "Chromosome", "Start_Position", "Reference_Allele",
#' and "Tumor_Seq_Allele2".
#' @param keep_vcf_id A \code{logical} value indicating whether to keep the
#' original ID field from the VCF file in the output. Default is \code{TRUE}.
#'
#' @return A \code{GRanges} object containing the variants loaded from the
#' specified VCF or MAF file. For VCF files, only the \code{ID}, \code{REF},
#' and \code{ALT} metadata columns are included in the output. For MAF files,
#' all available MAF columns are retained.
#'
#' @export
#'
#' @examples
#' # load the vcf file, please do not forget to specify the file type
#' vcf_path <- system.file("extdata",
#'                         "example_variants_chr16.vcf",
#'                         package = "G4SNVHunter")
#' variants <- loadVariant(vcf_path, file_type = "vcf")
#' # load the maf file
#' maf_path <- system.file("extdata",
#'                         "example_variants_chr16.maf",
#'                         package = "G4SNVHunter")
#' variants <- loadVariant(maf_path, file_type = "maf")
loadVariant <- function(variant_file,
                        file_type = c("vcf", "maf"),
                        keep_vcf_id = TRUE) {

  validateLoadVariantInputs(variant_file,
                            file_type,
                            keep_vcf_id)

  result <- switch(
    match.arg(file_type),
    vcf = loadVcf(variant_file, keep_vcf_id),
    maf = loadMaf(variant_file)
  )

  return(result)

}

#' validateLoadVariantInputs
#'
#' This function validates the input arguments for the \code{loadVariant}
#' function.
#'
#' @param variant_file A \code{character} string specifying the path to the
#' variant file, which can be in either VCF or MAF format.
#' @param file_type A \code{character} string specifying the type of the
#' input file: either \code{"vcf"} or \code{"maf}. This parameter must be
#' provided. If \code{"maf"} is specified, the following columns must be
#' present in the MAF file: "Chromosome", "Start_Position", "Reference_Allele",
#' and "Tumor_Seq_Allele2".
#' @param keep_vcf_id A \code{logical} value indicating whether to keep the
#' original ID field from the VCF file in the output. Default is \code{TRUE}.
#'
#' @return NULL
#'
#' @importFrom tools file_ext
#'
#' @noRd
validateLoadVariantInputs <- function(variant_file,
                                      file_type,
                                      keep_vcf_id) {

  if (length(file_type) != 1) {
    stop("'file_type' must be a single specified value: ",
         "either 'vcf' or 'maf'.")
  }

  if (!file_type %in% c("vcf", "maf")) {
    stop("Invalid file type: must be 'vcf' or 'maf'.")
  }

  file_ext <- tolower(file_ext(variant_file))
  if (file_ext == "vcf" && file_type == "maf") {
    stop("The file looks like a VCF (based on extension), but you ",
         "specified file_type = 'maf'.\nPlease verify and correct ",
         "either the file_type or the file extension.")
  } else if (file_ext == "maf" && file_type == "vcf") {
    stop("The file looks like a MAF (based on extension), but you ",
         "specified file_type = 'maf'.\nPlease verify and correct ",
         "either the file_type or the file extension.")
  }

  if (!file.exists(variant_file)) {
    stop("File does not exist: ", variant_file, ".")
  }

  if (!is.logical(keep_vcf_id) || length(keep_vcf_id) != 1) {
    stop("keep_vcf_id must be a single logical value.")
  }
}

#' loadVcf
#'
#' This function loads small variants such as SNVs, INDELs, and DELINs from a
#' VCF file. Structural variants and other types of variants will be ignored.
#'
#' @param vcf_file A \code{character} string specifying the path to the
#' VCF file.
#' @param keep_vcf_id A \code{logical} value indicating whether to keep the
#' VCF ID column in the output.
#'
#' @return A \code{GRanges} object containing the filtered small variant data.
#'
#' @importFrom VariantAnnotation readVcf isSNV isIndel isDelins
#' @importFrom VariantAnnotation isSubstitution isTransition
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom S4Vectors mcols
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @noRd
loadVcf <- function(vcf_file, keep_vcf_id) {

  vcf <- readVcf(vcf_file)
  initial_count <- length(vcf)
  if (initial_count == 0) stop("No variants found in the VCF file.")

  type_mat <- cbind(
    SNV = isSNV(vcf, singleAltOnly = FALSE),
    Indel = isIndel(vcf, singleAltOnly = FALSE),
    Delins = isDelins(vcf, singleAltOnly = FALSE),
    Substitution = isSubstitution(vcf, singleAltOnly = FALSE),
    Transition = isTransition(vcf, singleAltOnly = FALSE)
  )
  vcf <- vcf[rowSums(type_mat) > 0, ]
  filtered_count <- initial_count - length(vcf)
  if (length(vcf) == 0) stop("No variants passed the type filtering!")

  gr <- rowRanges(vcf)

  if (any(as.character(gr$REF) == "-") |
      any(unlist(as.character(unlist(gr$ALT))) == "-")) stop(
        "Detected MAF-style variants (with '-' notation) in the VCF file.")

  if (keep_vcf_id) {
    mcols(gr)$ID <- names(gr)
  }

  mcols(gr) <- mcols(gr)[, c(if(keep_vcf_id) "ID", "REF", "ALT")]

  message(sprintf(
    "VCF loading summary:
    - Initial variants: %d
    - Filtered out by type: %d
    - Final variants retained: %d",
    initial_count,
    filtered_count,
    length(gr)
  ))

  return(gr)
}

#' loadMaf
#'
#' This function loads small variants from a MAF file.
#' Structural variants and other types of variants will be ignored.
#'
#' @param maf_file A \code{character} string specifying the path to the
#' MAF file.
#'
#' @return A \code{GRanges} object containing the filtered small variant data.
#'
#' @importFrom data.table fread .SD
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @noRd
loadMaf <- function(maf_file) {

  maf <- fread(maf_file, header = TRUE, stringsAsFactors = FALSE)
  initial_count <- nrow(maf)

  required_cols <- c("Chromosome", "Start_Position",
                     "Reference_Allele", "Tumor_Seq_Allele2")
  missing_cols <- setdiff(required_cols, colnames(maf))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in MAF: ",
         paste(missing_cols, collapse = ", "))
  }

  maf <- maf[
    !is.na(get("Chromosome")) & get("Chromosome") != "" &
      !is.na(get("Start_Position")) & get("Start_Position") != "" &
      !is.na(get("Reference_Allele")) & get("Reference_Allele") != "" &
      !is.na(get("Tumor_Seq_Allele2")) & get("Tumor_Seq_Allele2") != "" &
      as.numeric(get("Start_Position")) == floor(
        as.numeric(get("Start_Position"))) &
      grepl("^[ACGTNacgtn-]+$", get("Reference_Allele")) &
      grepl("^[ACGTNacgtn-]+$", get("Tumor_Seq_Allele2")) &

      !startsWith(get("Reference_Allele"), "<") &
      !startsWith(get("Tumor_Seq_Allele2"), "<")
    , ]

  filtered_count <- initial_count - nrow(maf)

  message(sprintf(
    "MAF loading summary:
    - Initial variants: %d
    - Filtered out by type: %d
    - Final variants retained: %d",
    initial_count,
    filtered_count,
    nrow(maf)
  ))

  gr <- GRanges(
    seqnames = maf$Chromosome,
    ranges = IRanges(
      start = maf$Start_Position,
      end = maf$Start_Position + nchar(maf$Reference_Allele) - 1
    ),
    strand = "*",
    Reference_Allele = maf$Reference_Allele,
    Tumor_Seq_Allele2 = maf$Tumor_Seq_Allele2,
    maf[, .SD, .SDcols = setdiff(colnames(maf), required_cols)]
  )

  return(gr)
}

