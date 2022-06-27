
function  [yr,cntr,ynoiser]  =   parallel_3D_omp(I,params)

p = params.blocksize(3);

[k1,k2,N]  = size(I);
params.lambda = 0;
params.dict = [params.D1,params.D2];
params.select = size(params.dict,2);



x       = zeros(k1,k2,p,N-p+1,'single');
% y       = zeros(k1,k2,p,N-p+1,'single');
% ynoise  = zeros(k1,k2,p,N-p+1,'single');
cnt     = zeros(k1,k2,p,N-p+1,'single');

yr        =  zeros(size(I),'single');
cntr      =  zeros(size(I),'single');
% ynoiser   =  zeros(size(I),'single');
% parpool(4);
% tic
for i = 1:N-p+1
    x(:,:,:,i) = I(:,:,i:i+p-1);
end
clear I

parfor i = 1:N-p+1
% for i = 1:N-p+1
%     [y(:,:,:,i),cnt(:,:,:,i)] = newompdenoise3_4_parallel(params,x(:,:,:,i),5);
%       [x(:,:,:,i),cnt(:,:,:,i),ynoise(:,:,:,i)] = omp_3d_Distinctive_parallel(params,x(:,:,:,i),0);
      [x(:,:,:,i),cnt(:,:,:,i),~] = omp_3d_Distinctive_parallel(params,x(:,:,:,i),0);
%       fprintf( 'DFR-->>> slice: %2d\n', i);
end
% toc
% delete(gcp);

for i = 1:N-p+1
    yr(:,:,i:i+p-1)    = yr(:,:,i:i+p-1) + x(:,:,:,i);
    cntr(:,:,i:i+p-1)  = cntr(:,:,i:i+p-1) + cnt(:,:,:,i);
%     ynoiser(:,:,i:i+p-1)  = ynoiser(:,:,i:i+p-1) + ynoise(:,:,:,i);
   ynoiser   =1;
end




