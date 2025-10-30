# BiocManager::install("harmonicmeanp")
# BiocManager::install("EmpiricalBrownsMethod")

library(EmpiricalBrownsMethod)
library(harmonicmeanp)

# # ======= Example Data Just for Testing =======
# set.seed(123)
# # Create a mock 10x8 abundance matrix (10 kmers, 8 cells)
# kmer_abundance <- matrix(
#   rpois(80, lambda = 10),  # Poisson-distributed counts
#   nrow = 10,
#   ncol = 8,
#   dimnames = list(paste0("kmer_", 1:10), paste0("cell_", 1:8))
# )
# # Mock p-values for each kmer (from some prior statistical test)
# pvals <- runif(10, min = 0.0001, max = 0.2)
# ===============================

# ======= P-value Combination Methods =======
# Empirical Brown's Method
combine_pvalues_ebm <- function(pvals, data_matrix) {
  df <- as.data.frame(data_matrix)
  result <- EmpiricalBrownsMethod::empiricalBrownsMethod(
    data_matrix = df,
    p_values = pvals
  )
  return(result)
}

# Harmonic Mean P-value
combine_pvalues_hmp <- function(pvals) {
  result <- harmonicmeanp::p.hmp(pvals, L = length(pvals))
  return(result)
}

# ======= Example Test =======
cat("K-mer abundance matrix:\n")
print(kmer_abundance)

cat("\nMock p-values:\n")
print(pvals)

cat("\n--- Empirical Brown's Method ---\n")
ebm_result <- combine_pvalues_ebm(pvals, kmer_abundance)
print(ebm_result)

cat("\n--- Harmonic Mean P-value ---\n")
hmp_result <- combine_pvalues_hmp(pvals)
print(hmp_result)
