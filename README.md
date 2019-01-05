# multiscale-sampling-supplement

This repository contains Matlab code to accompany the article
> Krithika Manohar, Eurika Kaiser, Steven L. Brunton and J. Nathan Kutz. "Optimized Sampling for Multiscale Dynamics". *SIAM Multiscale Modeling and Simulation* (2019). To Appear.

A preprint of this article is available on [arXiv](https://arxiv.org/pdf/1712.05085.pdf). 

This work develops methods for optimal sampling in spatial domains (i.e., sensor placement) for discovering and estimating dynamics operating on multiple time scales. The primary tools are dimensionality reduction using multiresolution dynamic mode decomposition (mrDMD) and matrix QR pivoting for sensor placement.
In this code we provide the following:
- Functions to compute the mrDMD and optimal sensor placements for a given dataset.
- Scripts that generate the figures in this paper. 
- Examples of multiscale sampling on NOAA sea surface temperature (SST) data and an articial multiscale video example. 

The mrDMD method and subroutines in 'src/' are adapted from previous work
> J. Nathan Kutz, Steven L. Brunton, Bingni W. Brunton and Joshua L. Proctor. *Dynamic mode Decomposition: data-driven modeling of complex systems*. Vol. 149. SIAM (2016).

## External dependencies

- Matlab Signal Processing Toolbox
- [SPGL1](https://www.cs.ubc.ca/~mpf/spgl1/) solver for sparse systems (included)
- [export_fig](https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig) (included)

## Installation

Run setup.m to configure Matlab path, external packages and datasets.
