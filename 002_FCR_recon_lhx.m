


clear all;  
close all;

% addpath('D:\My_Matlab_code\LIUHANXI\Mfile');
% addpath('D:\My_Matlab_code\DFR_post\DFR');
% addpath('D:\My_Matlab_code\DFR_post\DFR\ksvdbox13');
% addpath('D:\My_Matlab_code\DFR_post\DFR\ksvdbox13\ompbox10');
% addpath('D:\My_Matlab_code\CT_DL\Dictionary_learning\Utilities');
cd ('C:\Program Files\MATLAB\R2021a\toolbox\ksvdbox13\private')

%%  Phantom
pMatrix_file = 'D:\UIH-CBCT-8-16\PMatrix_DD_0_SA_200_PN_601_BM_4_SL_2_CW.raw';
recon_folder = 'C:\Users\Liuhanxi\Desktop\3D\FCR';
sod = 750;sdd = 1250; dod = sdd -sod;
nu = 648;
nv = 474;
du = 0.616; dv = 0.616; 
nx = 512;ny = 512; nz = 387;
nz = 20;
dx = 0.449503;dy = dx; dz =dx;
rfov = nx*dx; zfov = nz*dz;
nview = 601;
first_Angle = -20;
angularArc = 200;
angles = first_Angle + (0:nview-1)*angularArc/(nview-1);
angles = angles./180*pi;
vol_geom = astra_create_vol_geom(ny,nx,nz,-rfov/2,rfov/2,-rfov/2,rfov/2,-zfov/2,zfov/2);
proj_geom = astra_create_proj_geom('cone',du,dv,nv,nu,angles,sod,dod);
proj_geom = geometry_correction(pMatrix_file,nu,proj_geom);
load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\low_dose_proj.mat')
% proj_geom.DetectorRowCount = 40;
ig          = image_geom('nx', nx,'ny', ny, 'nz', nz,'fov', 50);
%%
   load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\DFR_64_4.mat')
    ksvd_params.blocksize   =   [8, 8, 4];
    ksvd_params.memusage    =   'high';
    ksvd_params.stepsize    =   [1, 1, 1];
    ksvd_params.maxval      =   0.5;
    ksvd_params.maxatoms    =   8;
    ksvd_params.trainnum    =   1000;
    ksvd_params.dictsize    =   500;
    ksvd_params.parfor      =   1;
    ksvd_params.D1          =   dictA;
    ksvd_params.D2          =   dictB;
    clear dictA
% %% 
% %%plot the geometry
% astra_plot_geom(proj_geom);
% hold on
% astra_plot_geom(vol_geom,dx,'Magnification',1,'LineWidth',1,'Color','r');


%% ASTRA SIRT3D 
proj_id = astra_mex_data3d('create','-proj3d',proj_geom,proj);
rec_id = astra_mex_data3d('create','-vol',vol_geom);
cfg = astra_struct('SIRT3D_CUDA');
cfg.ProjectionDataId = proj_id;
cfg.ReconstructionDataId = rec_id;
alg_id = astra_mex_algorithm('create',cfg);
astra_mex_algorithm('iterate',alg_id,300);
rec_3D = astra_mex_data3d('get',rec_id);
astra_mex_data3d('delete', rec_id);
astra_mex_data3d('delete', proj_id);
%%
One3D            =  ones([nx,ny,nz],'single');
[proj_id, v1]    =  astra_create_sino3d_cuda(One3D, proj_geom, vol_geom);          
astra_mex_data3d('delete', proj_id);
[vol_id, Norimg]           =  astra_create_backprojection3d_cuda(v1, proj_geom, vol_geom); 
astra_mex_data3d('delete', vol_id);
%%
par.mask         =  ig.ones;
par.pixmax       =  0.5;
par.niter        =  20;  
par.proj_geom    =  proj_geom;
par.vol_geom     =  vol_geom;
par.Norimg       =  Norimg;
param.par        =  par;
param.max_it     =  50;
%%
alphan        =  0.0001:0.0001:0.0001;
lambdan       =  1000:200:1600;
res_path1     = 'C:\Users\Liuhanxi\Desktop\3D\FCR\';

for n = 1:length(lambdan)*length(alphan)    
    kk = mod(n-1, length(lambdan)) + 1;
    jj = floor((n-1)/length(lambdan)) + 1;   
param.mu              =  lambdan(kk);  
ksvd_params.sigma     =  alphan(jj);
param.ksvd            =  ksvd_params;
Recon_all  =  FCR_CT(rec_3D,proj,param);
NAME = ['FCR_Recon_',num2str(param.mu),'_',num2str(ksvd_params.sigma)];
save([res_path1,NAME,'.mat'], 'Recon_all');
end




