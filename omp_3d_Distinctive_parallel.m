



function [imout,weight,imnoise] = omp_3d_Distinctive_parallel(params,x,msgdelta)

%OMPDENOISE3 OMP for 3-D signals.
%
%  See also OMPDENOISE.

%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% parse input arguments %
if (nargin <2)
  x = params.x;
end


% D = params.dict;
D = [params.D1, params.D2];
blocksize = params.blocksize;


% blocksize %
if (numel(blocksize)==1)
  blocksize = ones(1,3)*blocksize;
end


% maxval %
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval = 1;
end


% gain %
if (isfield(params,'gain'))
  gain = params.gain;
else
  gain = 1.15;
end


% maxatoms %
if (isfield(params,'maxatoms'))
  maxatoms = params.maxatoms;
else
%   maxatoms = floor(prod(blocksize)/2);
end


% stepsize %
if (isfield(params,'stepsize'))
  stepsize = params.stepsize;
  if (numel(stepsize)==1)
    stepsize = ones(1,3)*stepsize;
  end
else
  stepsize = ones(1,3);
end
if (any(stepsize<1))
  error('Invalid step size.');
end


% noise mode %
if (isfield(params,'noisemode'))
  switch lower(params.noisemode)
    case 'psnr'
      sigma = maxval / 10^(params.psnr/20);
    case 'sigma'
      sigma = params.sigma;
    otherwise
      error('Invalid noise mode specified');
  end
elseif (isfield(params,'sigma'))
  sigma = params.sigma;
elseif (isfield(params,'psnr'))
  sigma = maxval / 10^(params.psnr/20);
else
  error('Noise strength not specified');
end


% lambda %
if (isfield(params,'lambda'))
  lambda = params.lambda;
else
  lambda = maxval/(10*sigma);
end


% msgdelta %
if (nargin <3)
  msgdelta = 5;
end
if (msgdelta<=0)
  msgdelta = -1;
end


epsilon = sqrt(prod(blocksize)) * sigma * gain;   % target error for omp


MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% compute G %

G = [];
if (memusage >= MEM_NORMAL)
  G = D'*D;
end


% verify dictionary normalization %

if (isempty(G))
  atomnorms = sum(D.*D);
else
  atomnorms = diag(G);
end
if (any(abs(atomnorms-1) > 1e-2))
  error('Dictionary columns must be normalized to unit length');
end

param.L=maxatoms;
param.eps= epsilon;
param.numThreads=-1;

nz = 0;  % count non-zeros in block representations

% the denoised signal
y = zeros(size(x),'single');
ynoise = zeros(size(x),'single');
cnt = zeros(size(x),'single');

blocknum = prod(floor((size(x)-blocksize)./stepsize) + 1);
processedblocks = 0;
tid = timerinit('ompdenoise', blocknum);

Jnumj = 20;
Jnumk = 10;
% Jnumj = 100;
% Jnumk = 1;
% tic;
% poolobj = gcp;
tt1 = 1;
for k = 1:stepsize(3):size(y,3)-blocksize(3)+1
  for j = 1:(Jnumj*stepsize(2)):size(y,2)-stepsize(2)-blocksize(2)+1

    jumpSizej = min(j+(Jnumj*stepsize(2))-1,size(y,2)-stepsize(2)-blocksize(2)+1);
    blocks    = im2colstep(double(x(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1)),blocksize,stepsize);


      [blocks, dc] = remove_dc(blocks,'columns');
      gamma = omp2(D'*blocks,sum(blocks.*blocks),G,epsilon,'maxatoms',maxatoms,'checkdict','off');
%       gamma = omp(D'*blocks,G,maxatoms,'checkdict','off');
%       gamma = mexOMP(blocks, D, param);
      
%       nz = nz + nnz(gamma);
         
         cn = gamma ;
         cn(1:size(params.D1,2),:) =0;
         blocks2 = D*cn;
         
%%
% tt = size(gamma,2);
% tt2  =  tt1+tt;
% % noisegamma(:,tt1:tt2-1) = gamma(size(params.D1,2)+1:end,:);
% noisegamma(:,tt1:tt2-1) = gamma;
% tt1  =  tt1+tt;



  %%       
%       [~, dc2] = remove_dc(blocks2,'columns');   
         
    gamma(size(params.D1,2)+1:end,:) = 0;
%       noise_gamma  = gamma(size(params.D1,2)+1:end,:);
%       ind = noise_gamma  > 0.03;
%       noise_gamma(ind) = 0;
%       gamma(size(params.D1,2)+1:end,:) = noise_gamma;

    cleanblocks = add_dc(D*gamma, dc, 'columns');
%     cleanblocks = add_dc(cleanblocks, dc2, 'columns');
%     cleanblocks = D*gamma;
%%
%    cleanblocks = add_dc((D.*params.activity)*gamma, dc, 'columns');
%    blocks2     = (D.*~params.activity)*gamma;
%%    
    cleanvol = col2imstep(cleanblocks,[size(y,1) blocksize(2)+jumpSizej-j blocksize(3)],blocksize,stepsize);
%     cleanvol = col2imstep(cleanblocks,[size(y,1) blocksize(2:3)],blocksize,stepsize);
%     y(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) = y(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) + cleanvol;
    y(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1) = y(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1) + cleanvol; 
    
%     cleanvol = col2imstep(blocks2,[size(y,1) blocksize(2:3)],blocksize,stepsize);
    cleanvol = col2imstep(blocks2,[size(y,1) blocksize(2)+jumpSizej-j blocksize(3)],blocksize,stepsize);
%     ynoise(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) = ynoise(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) + cleanvol;
    ynoise(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1) = ynoise(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1) + cleanvol; 
    
    cleanblocks_c = ones(size(blocks2));
    cleanim_c = col2imstep(cleanblocks_c,[size(y,1) blocksize(2)+jumpSizej-j blocksize(3)],blocksize,stepsize);
    cnt(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1) =  cnt(:,j:jumpSizej+1+blocksize(2)-2,k:k+blocksize(3)-1) + cleanim_c;  
  end
  
end


% delete(poolobj);
% toc;

if (msgdelta>0)
  timereta(tid, blocknum);
end


nz = nz/blocknum;  
imout = y ;
imnoise = ynoise;
weight = cnt;




