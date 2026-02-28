library(EmpiricalBrownsMethod)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# ======= Config =======
CONTIG_FILE     <- args[1]
KMER_PVAL_FILE  <- args[2]
KMER_COUNT_FILE <- args[3]
OUTPUT_FILE_EBM <- args[4]

KMER_SIZE <- 31

# ======= Helpers =======
generate_kmers <- function(sequence, k = KMER_SIZE) {
  seq_length <- nchar(sequence)
  if (seq_length < k) return(character(0))

  kmers <- sapply(1:(seq_length - k + 1), function(i) {
    substr(sequence, i, i + k - 1)
  })

  unique(kmers)
}
combine_pvalues_ebm <- function(pvals, data_matrix) {
  if (is.null(data_matrix) || nrow(data_matrix) == 0 || ncol(data_matrix) == 0) return(NA_real_)
  pvals <- as.numeric(pvals)
  if (any(!is.finite(pvals))) stop("Non-finite p-values supplied")
  if (any(pvals <= 0 | pvals >= 1)) stop("p-values must be in (0,1)")

  tryCatch(
    EmpiricalBrownsMethod::empiricalBrownsMethod(
      data_matrix = as.data.frame(data_matrix),
      p_values = pvals
    ),
    error = function(e) NA_real_
  )
}


contigs_df <- fread(CONTIG_FILE)

kmer_pval_df <- fread(KMER_PVAL_FILE)
kmer_lookup <- setNames(kmer_pval_df[["padj"]], kmer_pval_df$kmer)

kmer_count_df <- fread(KMER_COUNT_FILE)
count_cols <- setdiff(colnames(kmer_count_df), "kmer")

# ======= Process contigs =======
out <- vector("list", nrow(contigs_df))

for (i in seq_len(nrow(contigs_df))) {
  contig_seq <- contigs_df$contig[i]

  contig_kmers  <- generate_kmers(contig_seq, k = KMER_SIZE)
  matched_kmers <- intersect(contig_kmers, names(kmer_lookup))

  if (length(matched_kmers) == 0) {
    out[[i]] <- data.frame(contig_seq = contig_seq, pvalue = NA_real_)
    next
  }

  matched_pvals <- kmer_lookup[matched_kmers]

  idx <- match(matched_kmers, kmer_count_df$kmer)
  matched_counts <- kmer_count_df[idx, ..count_cols]

  ebm_pval <- combine_pvalues_ebm(matched_pvals, matched_counts)

  out[[i]] <- data.frame(contig_seq = contig_seq, pvalue = ebm_pval)
}

ebm_df <- rbindlist(out, use.names = TRUE, fill = TRUE)
fwrite(ebm_df, OUTPUT_FILE_EBM, sep = "\t", quote = FALSE)