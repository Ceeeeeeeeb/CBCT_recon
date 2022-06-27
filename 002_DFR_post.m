clear all;  
close all;
% addpath('D:\My_Matlab_code\DFR_post\DFR');
% addpath('D:\My_Matlab_code\DFR_post\DFR\ksvdbox13');
% addpath('D:\My_Matlab_code\DFR_post\DFR\ksvdbox13\ompbox10');
% addpath('D:\My_Matlab_code\CT_DL\Dictionary_learning\Utilities');
cd ('C:\Program Files\MATLAB\R2021a\toolbox\ksvdbox13\private')
%% load the dictionary 
load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\dictionary\DFR_64_4_SIRT.mat')
res_path1   = 'C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\reconstruction\';
res_name1   = 'DFR_post_884_SIRT.mat';
%% DFR parameters
    ksvd_params.blocksize   =   [8, 8, 4];
    ksvd_params.memusage    =   'high';
    ksvd_params.stepsize    =   [1, 1, 1];
    ksvd_params.maxval      =   10;
    ksvd_params.maxatoms    =   8;
    ksvd_params.trainnum    =   10000;
    ksvd_params.dictsize    =   500;
    ksvd_params.parfor      =   0;
    ksvd_params.D1          =   dictA(:,1:1:end);
    ksvd_params.D2          =   dictB(:,1:1:end);
    ksvd_params.sigma       =   0; 
    ksvd_params.lambda      =   0; 
%% load the projection or reconstruction 
load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\reconstruction\low_dose_SIRT_all.mat')
% X_fbp   =single(rec3D); 
X_fbp = single(rec_3D2);
%%
 parpool(8);   
 ksvd_params.x = single(X_fbp);
 [z1,~] = parallel_3D_ompBIG(ksvd_params.x,ksvd_params,60) ;
 delete(gcp);
     DFR  =  z1;
%%  save the result

save([res_path1,res_name1], 'DFR');

%%     
