# PopGeneticsVis
PopGeneticsVis: a toolkit for whole-genome resequencing and bisulfite sequencing downstream visualization in population genetics.

## 1. Installation

### 1.1 Download from GitHub Release

```bash
# 1. Test
Rscript ./PopGeneticsVis/bin/plink2_pca_vis.R --help

usage: ./PopGeneticsVis/bin/plink2_pca_vis.R [-h] --eigenvec EIGENVEC --eigenval EIGENVAL
                             --sample_pop SAMPLE_POP --select_pops SELECT_POPS
                             [--x_pc X_PC] [--y_pc Y_PC]
                             [--point_size POINT_SIZE]
                             [--point_alpha POINT_ALPHA]

PCA Visualization for PLINK PCA Eigenvec

optional arguments:
  -h, --help            show this help message and exit
  --eigenvec EIGENVEC   PLINK PCA results [pca.eigenvec].
  --eigenval EIGENVAL   PLINK PCA results [pca.eigenval].
  --sample_pop SAMPLE_POP
                        Table with samples in col1, pops in col2
                        [sample_pop.txt].
  --select_pops SELECT_POPS
                        Select populations ['GroupA,GroupB'].
  --x_pc X_PC           Which PC for X Axis [1].
  --y_pc Y_PC           Which PC for Y Axis [2].
  --point_size POINT_SIZE
                        Points size [4].
  --point_alpha POINT_ALPHA
                        Points alpha [0.8].
```

### 1.2 Clone Git Repository

```bash
# 1. Clone
git clone https://github.com/benben-miao/PopGeneticsVis.git

# 2. Test
Rscript ./PopGeneticsVis/bin/plink2_pca_vis.R --help

usage: ./PopGeneticsVis/bin/plink2_pca_vis.R [-h] --eigenvec EIGENVEC --eigenval EIGENVAL
                             --sample_pop SAMPLE_POP --select_pops SELECT_POPS
                             [--x_pc X_PC] [--y_pc Y_PC]
                             [--point_size POINT_SIZE]
                             [--point_alpha POINT_ALPHA]

PCA Visualization for PLINK PCA Eigenvec

optional arguments:
  -h, --help            show this help message and exit
  --eigenvec EIGENVEC   PLINK PCA results [pca.eigenvec].
  --eigenval EIGENVAL   PLINK PCA results [pca.eigenval].
  --sample_pop SAMPLE_POP
                        Table with samples in col1, pops in col2
                        [sample_pop.txt].
  --select_pops SELECT_POPS
                        Select populations ['GroupA,GroupB'].
  --x_pc X_PC           Which PC for X Axis [1].
  --y_pc Y_PC           Which PC for Y Axis [2].
  --point_size POINT_SIZE
                        Points size [4].
  --point_alpha POINT_ALPHA
                        Points alpha [0.8].
```
