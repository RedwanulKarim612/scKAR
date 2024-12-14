# scKAR

## Parameters
The method's parameters can be defined in the `src/config.env` file. The file will look like this:
```shell
INPUT_DIR="path/to/data"

#Parameters for clustering
CLUSTERING_ALGO=leiden
MIN_GENES=10
MIN_CELLS=3
N_NEIGHBORS=10
N_PCS=40
RESOLUTION=0.6
```
The `INPUT_DIR` parameter should have the absolute path to the dataset.
The dataset folder should have the following structure:
```
data/
├── reads/
│   ├── sample_1.fastq.gz
│   └── sample_2.fastq.gz
└── expression_matrix
    └── expression.csv
```

## How to run:
From the terminal run the following commands:
```
cd ./src
bash ./run.sh
```
