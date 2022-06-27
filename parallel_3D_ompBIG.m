

function    [Y ,Res]    =   parallel_3D_ompBIG(X,params,step)

c    =  size(X,3);
indx = 1:step:c-step;
indx = [indx,c+1];

k = params.blocksize(3);
IMouts   = ones(size(X),'single');
Weights  = ones(size(X),'single');
% Noise    = ones(size(X),'single');

    n           = 1;
    id1         = indx(n);
    id2         = indx(n+1)-1;
    tic;
    [IM,We,~]   = parallel_3D_omp(X(:,:,id1:id2),params);
    IMouts(:,:,1:id2-k+1)  =  IM(:,:,1:id2-k+1);
    Weights(:,:,1:id2-k+1) =  We(:,:,1:id2-k+1);
%     Noise(:,:,1:id2-k+1)   =  No(:,:,1:id2-k+1);
    timecost =  toc;    
    
fprintf( 'DL===>>> Iter:  1;  TIME: %3.4f\n',timecost);
%%
dd = 1;
for n = 2:length(indx)-2
    tic;
    id1 = indx(n)-k+1;
    id2 = indx(n+1)-1;
    
    [IM,We,~]               =  parallel_3D_omp(X(:,:,id1:id2),params);
    IMouts(:,:,id1:id2-k+1)  =  IM(:,:,1:end-k+1);
    Weights(:,:,id1:id2-k+1) =  We(:,:,1:end-k+1);
%     Noise(:,:,id1:id2-k+1)   =  No(:,:,1:end-k+1);
    dd = dd+1;
    
timecost =  toc;     
    fprintf( 'DL===>>> Iter: %2d;  TIME: %3.4f\n', n,timecost);
end
%%

     n = length(indx)-1;
 if n ~= 1;
    id1 = indx(n)-k+1;
    id2 = indx(n+1)-1;
    
    [IM,We,~] = parallel_3D_omp(X(:,:,id1:id2),params);    
    IMouts(:,:,id1:id2)  =  IM;
    Weights(:,:,id1:id2) =  We;
%     Noise(:,:,id1:id2)   =  No;
 end
     
  Y    =    (IMouts+params.lambda.*X)./(Weights+params.lambda);   Y(isnan(Y)) = 0;   
%   Res  =    (Noise + params.lambda.*X)./(Weights+params.lambda);  Res(isnan(Res)) = 0;   
    Res =0;
end




