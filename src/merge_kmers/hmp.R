library(harmonicmeanp)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

# ======= Configuration =======
CONTIG_FILE <- args[1]
KMER_FILE   <- args[2]
OUTPUT_FILE_HMP <- args[3]
KMER_SIZE <- 31

# ======= Helper Function =======
generate_kmers <- function(sequence, k = KMER_SIZE) {
  seq_length <- nchar(sequence)
  if (seq_length < k) return(character(0))

  kmers <- sapply(1:(seq_length - k + 1), function(i) {
    substr(sequence, i, i + k - 1)
  })

  unique(kmers)
}

# ======= Harmonic Mean P-value (HMP) =======
combine_pvalues_hmp <- function(pvals) {
  if (length(pvals) == 0) return(NA_real_)
  pvals <- as.numeric(pvals)

  eps_low <- 1e-200
  eps_high <- 1 - 1e-16 
  pvals[pvals < eps_low] <- eps_low
  pvals[pvals > eps_high] <- eps_high

  # sanity check
  if (any(!is.finite(pvals))) {
    stop("Non-finite p-values supplied to combine_pvalues_hmp")
  }

  tryCatch(
    harmonicmeanp::p.hmp(pvals, L = length(pvals)),
    error = function(e) {
      cat("Error in HMP calculation:", e$message, "\n")
      NA_real_
    }
  )
}

contigs_df <- fread(CONTIG_FILE, sep = "\t", header = TRUE)

cat("\nReading k-mer file:", KMER_FILE, "\n")
kmer_df <- fread(KMER_FILE, sep = "\t", header = TRUE)

# Create k-mer -> p-value lookup
kmer_lookup <- setNames(kmer_df[["padj"]], kmer_df$kmer)

# ======= Process Contigs =======
cat("\nProcessing contigs...\n")
results <- vector("list", nrow(contigs_df))

for (i in seq_len(nrow(contigs_df))) {
  contig_name <- contigs_df$name[i]
  contig_seq  <- contigs_df$contig[i]

  # Generate k-mers from contig
  contig_kmers <- generate_kmers(contig_seq, k = KMER_SIZE)

  # Find k-mers that exist in the k-mer file
  matched_kmers <- intersect(contig_kmers, names(kmer_lookup))

  if (length(matched_kmers) == 0) {
    results[[i]] <- data.frame(
      contig_seq  = contig_seq,
      hmp_pvalue  = 1
    )
    next
  }

  # Extract p-values for matched k-mers
  matched_pvals <- kmer_lookup[matched_kmers]
  hmp_pval <- combine_pvalues_hmp(matched_pvals)

  results[[i]] <- data.frame(
    contig_seq  = contig_seq,
    hmp_pvalue  = hmp_pval
  )
}

results_df <- do.call(rbind, results)

# ======= Save Results =======
hmp_df <- data.frame(
  contig = results_df$contig_seq,
  pvalue = results_df$hmp_pvalue
)

cat("\nSaving HMP results to:", OUTPUT_FILE_HMP, "\n")
fwrite(hmp_df, OUTPUT_FILE_HMP, sep = "\t", quote = FALSE)

cat("\nResults saved successfully!\n")