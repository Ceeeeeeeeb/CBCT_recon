%%Initialize
clear;
close all;
clc;

%%hyper param
foldername = 'D:\My_Matlab_code\LIUHANXI\low_dose\tif\';
pMatrix_file = 'C:\Users\Liuhanxi\Desktop\Recon_matlab\PMatrix_DD_0_SA_200_PN_601_BM_4_SL_2_CW.raw';
recon_folder = 'D:\My_Matlab_code\LIUHANXI\Results\Iter\';
% [0~1] double, 0 means no truncation correction
%truncation_weight = 1.0;
% filter_type = 'hann';
% filter_cutoff = 0.9;

%%define Geometry
sod = 750;sdd = 1250; dod = sdd -sod;
nu = 648;
nv = 474;


du = 0.616; dv = 0.616; 
nx = 512;ny = 512; nz = 387;
% nz = 20;
dx = 0.449503;dy = dx; dz =dx;
rfov = nx*dx; zfov = nz*dz;
nview = 601;
first_Angle = -20;
angularArc = 200;
angles = first_Angle + (0:nview-1)*angularArc/(nview-1);
angles = angles./180*pi;


iteration = 500;
% load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\low_dose_proj.mat')
load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\projection\high_dose_proj.mat')
% proj = DFR;
% proj_file_list = dir([foldername,'*.tif']);
% proj = zeros(648, nview, 474, 'single');
% for i=1:nview
%     proj(:,i,:) = imread([foldername, proj_file_list(i).name])';
% end    
%

%% astra geom & geometry correction
vol_geom = astra_create_vol_geom(ny,nx,nz,-rfov/2,rfov/2,-rfov/2,rfov/2,-zfov/2,zfov/2);
proj_geom = astra_create_proj_geom('cone',du,dv,nv,nu,angles,sod,dod);
proj_geom = geometry_correction(pMatrix_file,nu,proj_geom);


% %%plot the geometry
% astra_plot_geom(proj_geom);
% hold on
% astra_plot_geom(vol_geom,dx,'Magnification',1,'LineWidth',1,'Color','r');

%recon

% proj =proj(:,:,218:257);
% proj_geom.DetectorRowCount = 40;
%% ASTRA SIRT3D 
proj_id = astra_mex_data3d('create','-proj3d',proj_geom,proj);
rec_id = astra_mex_data3d('create','-vol',vol_geom);
cfg = astra_struct('SIRT3D_CUDA');
cfg.ProjectionDataId = proj_id;
cfg.ReconstructionDataId = rec_id;
alg_id = astra_mex_algorithm('create',cfg);
astra_mex_algorithm('iterate',alg_id,iteration);
rec_3D = astra_mex_data3d('get',rec_id);
astra_mex_data3d('delete', rec_id);

%FOV crop
cropR = (nu*du/2*sod)/sdd;
maxD = min(nx-1)/2;
cropR = min([cropR/dx maxD]);
[x,y] = meshgrid(1:nx,1:ny);
inM = (x-nx/2).^2 + (y-ny/2).^2<cropR^2;
rec_3D2 = bsxfun(@times, rec_3D ,inM);


%% Save the result
% wl = 0.02;
% ww = 0.02;
% figure, imshow(fliplr(squeeze(rec_3D(:,:,194))),[wl-ww/2, wl+ww/2]);
% img = single(squeeze(rec_3D(:,:,194)))';
% filename = fullfile(recon_folder,'recon_194.tif');
% fTIF = Fast_Tiff_Write(filename);
% fTIF.WriteIMG(img);
% fTIF.close;
save('high_dose_SIRT_all','rec_3D2')
% for i = 1:nz
%     img = single(fliplr(squeeze(rec_3D(:,:,i))))';
%     fTIF =Fast_Tiff_Write([recon_folder, 'slice_', num2str(i, '%04d'), '.tif']);
%     fTIF.WriteIMG(img);
%     fTIF.close;
% end


%% Clean up
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', rec_id);
astra_mex_data3d('delete', proj_id);






