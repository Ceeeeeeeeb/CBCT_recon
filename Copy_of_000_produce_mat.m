clear
close all;

%% saving the projection or reconstruction
is_proj = 0;

%%
foldername = 'C:/Users/Liuhanxi/Desktop/Recon_matlab/CBCT_Scatter_Correction-main/CBCT_recon/example/MRCPs/recon/phantom 3/x1/proj_primary_true/';
proj_file_list = dir([foldername,'*.tif']);

if is_proj==1
    nview = 401;
    nu = 324;
    nv = 237;
    proj = zeros(nu,nview,nv,'single');
    for i = 1:nview
        proj(:,i,:) = imread([foldername,proj_file_list(i).name])';
    end
    save('low_dose_proj_monte','proj')
else
    nview = 387;
    nu = 512;
    nv = 512;
    recon = zeros(nu,nview,nv,'single');
    for i = 1:nview
        recon(:,i,:) = imread([foldername,proj_file_list(i).name])';
    end
    save('high_dose_recon_monte','recon')
end

