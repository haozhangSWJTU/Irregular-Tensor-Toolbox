# Irregular-Tensor-Toolbox (Matlab)
Hao Zhang, Ting-Zhu Huang, Xi-Le Zhao, Shuqin Zhang, Jin-Yu Xie, Tai-Xiang Jiang, Michael K. Ng

<p align="center">
    <br>
    <img src="flowchart.png" width="60%"/>
    <br>
</p>

## Abstract
Tensor decompositions have been successfully applied to multidimensional data recovery. However, classical tensor decompositions are not suitable for emerging spatio-irregular multidimensional data (i.e., spatio-irregular tensor), whose spatial domain is non-rectangular, e.g., spatial transcriptomics data from bioinformatics and semantic units from computer vision. By using preprocessing (e.g., zero-padding or element-wise 0-1 weighting), the spatio-irregular tensor can be converted to a spatio-regular tensor and then classical tensor decompositions can be applied, but this strategy inevitably introduces bias information, leading to artifacts. How to design a tensor-based method suitable for emerging spatio-irregular tensors is an imperative challenge. To address this challenge, we propose a learnable transform-assisted tensor singular value decomposition (LTA-TSVD) for spatio-irregular tensor recovery, which allows us to leverage the intrinsic structure behind the spatio-irregular tensor. Specifically, we design a learnable transform to project the original spatio-irregular tensor into its latent spatio-regular tensor, and then the latent low-rank structure is captured by classical TSVD on the resulting regular tensor. Empowered by LTA-TSVD, we develop spatio-irregular low-rank tensor completion (SIR-LRTC) and spatio-irregular tensor robust principal component analysis (SIR-TRPCA) models for the spatio-irregular tensor imputation and denoising respectively, and we design corresponding solving algorithms with theoretical convergence. Extensive experiments including the spatial transcriptomics data imputation and hyperspectral image denoising show SIR-LRTC and SIR-TRPCA are superior performance to competing approaches and benefit downstream applications.

## Citation
@ARTICLE{SIRTD,
author = {Zhang, Hao and Huang, Ting-Zhu and Zhao, Xi-Le and Zhang, Shuqin and Xie, Jin-Yu and Jiang, Tai-Xiang and Ng, Michael K.},
title = {Learnable Transform-Assisted Tensor Decomposition for Spatio-Irregular Multidimensional Data Recovery},
journal = {ACM Trans. Knowl. Discov. Data},
year = {2024},
doi = {10.1145/3701235}

