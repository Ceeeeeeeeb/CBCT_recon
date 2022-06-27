# Dictionary learning for CBCT low dose denoise in projection domain or image domain
## to use this code
### 0. saving the orignal files of projections or reconstructions as MATLAB binary files using *produce_mat.m*
### 1. loading the low/high dose projection/reconstruction binary file pair and generate dictionary using *Dicts_construct.m*
### 2. process the low dose projection/reconstruction with the generated dictionary
> directly process the projection/reconstruction with OMP algorithum using *DFR_post.m*
> 
> iteratively reconstruc the processed projections with FCR algorithum using *FCR_recon.m* 
