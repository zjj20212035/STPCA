# STPCA-MP
Author: Junjing Zheng

Email: zjj20212035@163.com

> **Paper**: [Orientation-Aware Sparse Tensor PCA for Efficient Unsupervised Feature Selection](arXiv Link) | **Authors**: [Junjing Zheng, Xinyu Zhang, Weidong Jiang, Xiangfeng Qiu, Mingjian Ren]

 This code is for reproduction of Table 7 (Metrics of the comparative methods on real-world dataset) and Figure 7 (Visualization on PIE, JAFFE, and BreastMNIST) in our paper, which conducts unsupervised feature selection (UFS) and clustering.
If you are going to using our code, please cite our paper with the following:

## Project Structure
```tree
Project/
├── code/        # code for reproduction
|   ├──Evaluation  # evaluation metrics computation
|   ├──methods     # UFS algorithm
|   ├──utils       # tools for efficient experiment
|   ├──config_Figure7.m   # Configuration for Figure7
|   ├──config_Table7.m   # Configuration for Table7
|   ├──exp_Figure7_visualization # Visualization for Figure 7
|   ├──exp_Table7_clustering_POC # Main experiment for Table 7
├── data/               # datasets(.mat)   
└── README.md           # This file
```
## Environment and Setup
- **Operating System**: Windows 11 (may work on other Windows versions)
- **MATLAB**: Version R2023b or higher
- **Required Toolbox**: 
  - Parallel Computing Toolbox (for parallel processing)
  - Tensor Toolbox v3.5 (included in `utils/` folder)
- **Environment Preparation**
```matlab
% include tool path
addpath(genpath([pwd,'\methods']))
addpath(genpath([pwd,'\Evaluation']))
addpath(genpath([pwd,'\utils']))
% start parallel computation toolbox
try 
    parpool;
catch
    disp('parpool already available')
end
```
The above commands are also included in 'exp_Figure7_visualization' and 'exp_Table7_clustering_POC'.

## Data Preprocessing
- All datasets are included in the `data/` folder in `.mat` format.
- Each dataset contains the following common attributes:
  - `data` $d_1 \times d_2 \times n$ tensor ($n$=samples, $d_i$=$i$-mode dimensionality)
  - `label` $n\times 1$ class vector
  - `tensor_based` (0: non-tensor data, 1: tensor data)
  - `tensor_type` (either "tube-wise" or "slice-wise")
  - `tensor_size` (dimensionality of the feature space)
- No manual data preprocessing needed - all preprocessing is handled automatically by:
  - `config_Table7.m` (for Table 7 experiments)
  - `config_Figure7.m` (for Figure 7 experiments)

in which `normalization` is carried out by dividing each element by the maximum absolute value, scaling all values to the interval $[-1,1]$.

## Quick Reproduction
Open this project in MATLAB R2023b and
- **To reproduce Table 7 results**: 
```matlab
cd code
run('exp_Table7_clustering_POC.m');
```
- **To reproduce Figure 7 results**: 
```matlab
cd code
run('exp_Figure7_visualization.m');
```
## Parameters Introduction
- Parameters in STPCA-MP:
  - $\lambda$: regularization parameter of $\ell_{2,1}$-norm
  - $\eta$: regularization parameter of trace function
  - $\mathbf{M}$: specify an invertible matrix for conducting $\star_{\mathbf{M}}$
  - `orientation`: specify an order set for data orientation
- Parameters for experiments:
  - `grid_search`: the range of grid search for $\lambda$ and $\eta$
  - `N`: the repeated times of K-means
  - `NumFS`: the number of selected features

Due to K-means's sensitivity to intialization, the reproduction results could be slightly different to the results in our paper but within a reasonable range.

## Version
```bash
git clone https://github.com/zjj20212035/STPCA.git
cd STPCA
git checkout v1.1
```