#' Convert HDL LD Reference Panels to Text Format
#'
#' Converts the R-serialized (.rda) HDL reference panels to a plain-text
#' directory format that the gsem Rust CLI can read. This is a one-time
#' conversion step needed before running \code{genomicsem hdl}.
#'
#' @param ld.path Path to the HDL reference panel directory containing
#'   \code{UKB_snp_counter_*.rda}, \code{UKB_snp_list_*.rda}, and
#'   per-piece \code{chr*.rda} and \code{chr*.bim} files.
#' @param out.path Path to the output directory for text-format files.
#'   Created if it does not exist.
#' @param verbose Print progress messages (default TRUE).
#' @return Invisibly returns the output directory path.
#'
#' @details
#' The output directory will contain:
#' \describe{
#'   \item{pieces.tsv}{Tab-delimited index: chr, piece, n_snps}
#'   \item{chr\{N\}.\{P\}.snps.tsv}{Per-piece SNP data: SNP, A1, A2, LDsc}
#' }
#'
#' These files are read by \code{genomicsem hdl --ld-path <out.path>}.
#'
#' @examples
#' \dontrun{
#' # Download HDL reference panels (see HDL documentation)
#' # Then convert once:
#' convert_hdl_panels(
#'   ld.path = "path/to/UKB_LD_reference",
#'   out.path = "path/to/hdl_text_panels"
#' )
#'
#' # Now use with gsem CLI:
#' # genomicsem hdl --traits t1.gz t2.gz --ld-path hdl_text_panels/ -o hdl.json
#' }
#' @export
convert_hdl_panels <- function(ld.path, out.path, verbose = TRUE) {
  if (!dir.exists(ld.path)) {
    stop("LD reference directory not found: ", ld.path)
  }
  dir.create(out.path, showWarnings = FALSE, recursive = TRUE)

  # Load SNP counter to discover pieces
  counter.file <- list.files(ld.path, pattern = "UKB_snp_counter.*\\.(rda|RData)$",
                             full.names = TRUE)
  if (length(counter.file) == 0) {
    stop("UKB_snp_counter file not found in ", ld.path)
  }

  env <- new.env()
  load(counter.file[1], envir = env)
  nsnps.list <- env$nsnps.list.imputed

  # Build pieces index
  pieces <- data.frame(chr = integer(), piece = integer(), n_snps = integer())
  for (chr in seq_along(nsnps.list)) {
    for (p in seq_along(nsnps.list[[chr]])) {
      pieces <- rbind(pieces, data.frame(
        chr = chr, piece = p, n_snps = nsnps.list[[chr]][p]
      ))
    }
  }

  write.table(pieces, file.path(out.path, "pieces.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  if (verbose) message("Found ", nrow(pieces), " LD pieces across ",
                       length(nsnps.list), " chromosomes")

  # Convert each piece
  n_converted <- 0
  for (i in seq_len(nrow(pieces))) {
    chr <- pieces$chr[i]
    p <- pieces$piece[i]

    # Find the .rda file for this piece
    rda.pattern <- sprintf("chr%d\\.%d\\.", chr, p)
    rda.files <- list.files(ld.path, pattern = paste0(rda.pattern, ".*\\.(rda|RData)$"),
                            full.names = TRUE)
    bim.files <- list.files(ld.path, pattern = paste0(rda.pattern, ".*\\.bim$"),
                            full.names = TRUE)

    if (length(rda.files) == 0 || length(bim.files) == 0) {
      if (verbose) message("  Skipping chr", chr, " piece ", p, " (files not found)")
      next
    }

    # Load LD scores
    piece.env <- new.env()
    load(rda.files[1], envir = piece.env)
    ldsc <- piece.env$LDsc

    # Load bim file (PLINK format: CHR, SNP, CM, BP, A1, A2)
    bim <- read.table(bim.files[1], header = FALSE, stringsAsFactors = FALSE)
    colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

    n_snps <- min(length(ldsc), nrow(bim))
    out.df <- data.frame(
      SNP = bim$SNP[1:n_snps],
      A1 = bim$A1[1:n_snps],
      A2 = bim$A2[1:n_snps],
      LDsc = ldsc[1:n_snps]
    )

    out.file <- file.path(out.path, sprintf("chr%d.%d.snps.tsv", chr, p))
    write.table(out.df, out.file, sep = "\t", row.names = FALSE, quote = FALSE)
    n_converted <- n_converted + 1

    if (verbose && n_converted %% 50 == 0) {
      message("  Converted ", n_converted, " / ", nrow(pieces), " pieces")
    }
  }

  if (verbose) message("Done. Converted ", n_converted, " pieces to ", out.path)
  invisible(out.path)
}
