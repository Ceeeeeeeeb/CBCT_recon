clear
close all;

%% saving the projection or reconstruction
is_proj = 0;

if is_proj==1
    nview = 401;
    nu = 324;
    nv = 237;
else
    nview = 387;
    nu = 512;
    nv = 512;
end

foldername = 'C:/Users/Liuhanxi/Desktop/Recon_matlab/CBCT_Scatter_Correction-main/CBCT_recon/example/MRCPs/recon/phantom 3/x0.25/proj_primary_true/';
proj_file_list = dir([foldername,'*.tif']);
proj = zeros(nu,nview,nv,'single');

for i = 1:nview
    proj(:,i,:) = imread([foldername,proj_file_list(i).name])';
end

save('low_dose_proj_monte','proj')