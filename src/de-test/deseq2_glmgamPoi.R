suppressPackageStartupMessages({
    library(DESeq2)
    library(BiocParallel)
    library(ashr)
    library(data.table)
    library(scran)
    library(ggplot2)
})


register(MulticoreParam(n_cpus)) 

args <- commandArgs(trailingOnly = TRUE)

metadata_file <- args[1]
matrix_file <- args[2]
tpm_file <- args[3]
save_path <- args[4]
rowThreshold <- as.numeric(args[5])
colThreshold <- as.numeric(args[6])
n_cpus <- as.numeric(args[7])

tpm_data <- read.csv(tpm_file, header = FALSE, row.names = 1, sep = "\t")
colnames(tpm_data) <- c("tpm", "cluster")
total_tpm_sums <- tpm_data$tpm

sample_info <- read.csv(metadata_file, header = FALSE, sep = "\t")
colnames(sample_info) <- c("file", "condition")
sample_info$condition <- factor(sample_info$condition)

print(table(sample_info$condition))

# sample_info$condition <- droplevels(sample_info$condition)

file <- matrix_file

    print("--------------------------------------------")
    print("reading count matrix")
    count_data <- read.csv(file, header = FALSE, row.names = 1, fill = TRUE, na.strings = "")
    count_data[is.na(count_data)] <- 0
    colnames(count_data) <- sample_info$file

    # Unnormalizing TPM data
    for (i in 1:length(total_tpm_sums)) {
        count_data[,i] <- as.integer(round((count_data[,i] * total_tpm_sums[i]) / 1e6, digits = 0))
    }
    count_data[is.na(count_data)] <- 0

    # Filter out lowly expressed genes/kmers
    print("Filtering out lowly expressed genes/kmers")
    print(paste0("Number of rows before filtering: ", nrow(count_data)))
    keep <- rowSums(count_data)  # Adjust thresholds as needed
    # print(head(keep >= 100))
    # print(head(keep))
    count_data <- count_data[keep >= rowThreshold, ]
    print(paste0("Number of rows after filtering: ", nrow(count_data)))

    print(paste0("Number of columns before filtering: ", ncol(count_data)))
    keep <- colSums(count_data)
    # print(head(keep >= 5))
    count_data <- count_data[ ,keep >= colThreshold]
    sample_info <- sample_info[keep >= colThreshold, ]
    
    print(paste0("Number of columns after filtering: ", ncol(count_data)))

    rownames(sample_info) <- sample_info$file

    print(paste0("Starting DESeq2 analysis for file"))
    dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~ condition)

    rm(count_data); gc()

    print("Computing Sum factors (estimateSizeFactors)")
    dds <- computeSumFactors(dds, min.mean = 1e-6, cluster = sample_info$condition, positive = TRUE, BPPARAM = MulticoreParam(n_cpus))

    dds$condition <- factor(dds$condition, levels = c("A","B"))
    dds$condition <- droplevels(dds$condition)

    print("Running Deseq2 analysis")
    dds <- DESeq(dds, parallel = FALSE, test = "LRT", fitType = "glmGamPoi", minmu=1e-3, minRep=Inf, reduced = ~ 1)

    print(paste0('finish DESeq2 analysis for file: ', matrix_file))
    
    deseq_results <-  results(dds, name = "conditionB", parallel = TRUE)

    print(paste0('start writing results for file: ', matrix_file))

    write.table(as.data.frame(deseq_results), 
            file = save_path,
            sep = "\t", 
            row.names = TRUE, 
            col.names = !file.exists(save_path))
   
    print(paste0('filtered results for file:' , matrix_file, " has ", nrow(deseq_results), " rows"))
# }
# --------------------------------------------------------------------------------------------------------
