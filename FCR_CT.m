



function [Iout1,Iout2, param1] = FCR_CT(Image_int,Proj,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu        = param.mu;
num       = param.max_it;
par       = param.par;
x      = Image_int;
z      = Image_int;

ksvd_params = param.ksvd;

%%
reset(gpuDevice);
for kp=1:num    
tic;  
%% subproblem-1-2 %%%%%%%%%%%%%%%%%%
    x         =  SART_Cone_dict(x,par.proj_geom,par.vol_geom,Proj,par.Norimg,z,ones(size(x)),mu,par);
%% subproblem-3 %%%%%%%%%%%%%%%%%%    
                      [IMout,Weight,~] =  parallel_3D_omp(x,ksvd_params);
                      z = IMout./Weight; 
                      z(isnan(z)) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %Compute various performance metrics
%     param1.SSIM(kp) = Ssim(double(x),param.recon_ref);
%     param1.PSNR(kp) = psnr(double(x),param.recon_ref);
%     param1.RMSE(kp) = RMSE(double(x),param.recon_ref)/0.38*2000;
    time  = toc; 
%   fprintf( 'DL-CT-------->>>>>>Iter: %3d;  SSIM: %3.4f;PSNR: %3.4f;RMSE: %3.4f;  Time: %3.2f\n', kp, param1.SSIM(kp),param1.PSNR(kp),param1.RMSE(kp),time);
    fprintf( 'DL-CT-------->>>>>>Iter: %3d; Time: %3.2f\n', kp, time);  

end
Iout1 = x;
Iout2 = z;
end
