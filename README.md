# scKAR: Single-Cell K-mer Based Analysis and Reconstruction

scKAR is a reference free bioinformatics tool designed for k-mer-based differential expression analysis single-cell RNA sequencing (scRNA-seq) data. It allows detection of differential expression transcripts. The tool is focused for non-model organisms where the reference genome is scant and diseases where atypical splicing causes various genetic variation (intron retention, lncRNA/miRNA driven regulation), causing isoform expression which are not annotated in reference transcriptome. 

## Table of Contents
- [Installation](#installation)
- [Input Data Format](#input-data-format)
- [Configuration File](#configuration-file)
- [Reference Transcriptome Processing](#reference-transcriptome-processing)
- [Running scKAR](#running-sckar)
- [Output Files](#output-files)
- [Operation Modes](#operation-modes)
- [Mode of Condition Generation](#mode-of-condition-generation)
  
---

## Installation

scKAR requires Python and R dependencies. Standard C++ compiler is also needed. To install the required environment, follow these steps:

### Setting up Python Environment
```sh
conda env create -f python_env.yml
```

### Setting up R Environment
```r
install.packages("renv")
renv::restore()
```

---

## Input Data Format

The input dataset should be structured as follows:

```
data/
├── reads/
│   ├── sample_1.fastq.gz
│   └── sample_2.fastq.gz
├── expression_matrix/
|   └── expression.csv
├── metadata/
|   ├── clustering.csv
|   └── bipartitions.csv
└── genome_assembly/
    └── genome.2bit

```

- `reads/` directory contains the raw sequencing reads in FASTQ format.
- `expression_matrix/` contains the gene expression matrix in CSV format.
- `metadata/` contains necessary files to generate conditions when `MODE="CUSTOM_METADATA"`.
- `genome_assembly/` contains the complete genome assembly for the organism being studied.

The absolute path to this dataset should be provided in the `config.env` file.

---

## Configuration File

The `config.env` file contains all hyperparameters required to run scKAR. Below is an explanation of each parameter:

### Paths for Input Data
```sh
INPUT_DIR="/path/to/data/"
REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH="/path/to/data"
NUMBER_OF_THREADS=8
```
- `INPUT_DIR`: Absolute path to the dataset directory.
- `REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH`: Path to the preprocessed 31-mers reference transcriptome (required if `FILTER_REFERENCE_KMERS=TRUE`).
- `NUMBER_OF_THREADS`: 8, the number of parallel cores to use for filtering, adjacency matrix creation, F-test and DE-test.
  
### K-mer Space
```sh
FILTER_REFERENCE_KMERS=TRUE
```
- `FILTER_REFERENCE_KMERS`: If `TRUE`, the tool removes reference k-mers before analysis.

### Mode of Condition Generation
```sh
MODE="GENE_EXPRESSION"
# MODE="KMER_ABUNDANCE"
# MODE="CUSTOM_METADATA"
```
- `GENE_EXPRESSION`: Clustering based on gene expression.
- `KMER_ABUNDANCE`: Clustering based on k-mer abundance.
- `CUSTOM_METADATA`: Uses a predefined clustering and metadata.

### Parameters for Clustering
```sh
CLUSTERING_ALGO=leiden
MIN_GENES=10
MIN_CELLS=3
N_NEIGHBORS=10
N_PCS=40
RESOLUTION=0.6
```
- `CLUSTERING_ALGO`: Clustering algorithm (default: Leiden).
- `MIN_GENES`: Minimum number of genes expressed in a cell.
- `MIN_CELLS`: Minimum number of cells expressing a gene.
- `N_NEIGHBORS`: Number of nearest neighbors for graph construction.
- `N_PCS`: Number of principal components for dimensionality reduction.
- `RESOLUTION`: Resolution parameter for clustering.

### Parameters for F-Test
```sh
MIN_ROW_THREDSHOLD=40
PSEUDOBULK_SUZE=500000
```
- `MIN_ROW_THREDSHOLD`: Minimum threshold for row filtering before applying Fisher’s F-test.
- `PSEUDOBULK_SIZE`: The number of kmers that are pulled together for parallelized F-test. It is used for resource constraint mainly. If the matrix is too large, then independent pseudobulking is applied for augmenting the matrix into smaller matrices for processing. 

### Parameters for Differential Expression
```sh
MIN_ROW_COUNT=10
MIN_COL_COUNT=10
```
- `MIN_ROW_COUNT`: Minimum count of rows to consider a feature for differential expression.
- `MIN_COL_COUNT`: Minimum count of columns for valid feature consideration.

### Parameters for Final Assembly
```sh
LOG2FC=1.0
PVAL=0.05
BASE_MEAN_THRESHOLD=0.0
```
- `LOG2FC`: Log2 fold-change threshold for differential expression.
- `PVAL`: p-value threshold for significance testing.
- `BASE_MEAN_THRESHOLD`: Minimum mean expression threshold.

---

## Reference Transcriptome Processing

Before running scKAR, a preprocessed list of 31-mers from the reference transcriptome is required if `FILTER_REFERENCE_KMERS=TRUE`. To generate this list, execute:
```sh
bash preprocessing_script/sorted_31mers_generator_RT.sh
```
This generates a sorted list of 31-mers, which should be provided in `REFERENCE_TRANSCRIPTOME_SORTED_31MERS_PATH` in `config.env`.

---

## Running scKAR

To execute scKAR, run the following commands from the terminal:
```sh
cd ./src
bash ./run.sh
```

---

## Output Files

scKAR generates multiple output files categorized as follows:

### `bipartitions/`
- Stores metadata for differential expression (DE) testing, where each file contains cell-wise conditions ('A' or 'B').

### `deseq_results/`
- Intermediate DE test results for each pseudo-bulk matrix.

### `f_test_results/`
- Fisher’s F-test results stored in an adjacency matrix format.

### `final_results/`
Contains subdirectories for each bipartition, including:
```
bipartition_0/
├── A_kmers.fasta  # K-mers upregulated in condition A
├── B_kmers.fasta  # K-mers upregulated in condition B
└── abyss/
    ├── A_contigs.fasta  # Assembled contigs for condition A
    ├── B_contigs.fasta  # Assembled contigs for condition B
```
These files provide the final assembled contigs for downstream analysis.

---

## Operation Modes

scKAR can operate in different modes by modifying the `config.env` file:

- `FILTER_REFERENCE_KMERS=TRUE`: Removes reference k-mers before analysis.
- `MODE="GENE_EXPRESSION"`: Clustering based on gene expression.
- `MODE="KMER_ABUNDANCE"`: Clustering based on k-mer abundance.
- `MODE="CUSTOM_METADATA"`: Uses custom clustering and metadata.

---

## Mode of Condition Generation

If `MODE="CUSTOM_METADATA"`, users must provide the following files in `/data/metadata/` folder:

### `clustering.csv`
```
GSM1887215,1
GSM1887216,1
GSM1887217,0
GSM1887218,2
```
Format: `<cell_annotation>,<cluster_id>` (barcode/read file name, cluster ID).

### `bipartitions.csv`
```
set1	set2
{1, 2}	{0}
{1, 0}  {2}
```
Defines biconditions for differential expression analysis.

---

This README provides all necessary details to configure and run scKAR efficiently. If you encounter any issues, please raise an issue in the repository.

