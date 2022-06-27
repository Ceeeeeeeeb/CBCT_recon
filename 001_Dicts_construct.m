clear

% load('C:\Users\Liuhanxi\Desktop\Recon_matlab\FDK_lowdose.mat')
load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\reconstruction\monte\low_dose_recon_monte.mat')
%Rec_SART  =single(rec3D);
low_dose = single(recon);

% load('C:\Users\Liuhanxi\Desktop\Recon_matlab\FDK_highdose.mat')
load('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\reconstruction\monte\high_dose_recon_monte.mat')
% xtrue  =single(rec3D);
high_dose = single(recon);
Rec_SART = low_dose;
xtrue = high_dose;
%%
%512*512*387
%%
% xtrue      = double(xtrue(:,:,201:250));
xtrue      = double(xtrue(:,:,:));
Rec_SART   = double(Rec_SART(:,:,:));

Noise  =  xtrue -Rec_SART;
par.b    = [8,8,4];
par.slid = [2,2,2];

cd ('C:\Program Files\MATLAB\R2021a\toolbox\ksvdbox13\private')

DictsA= single(im2colstep(double(xtrue), par.b, par.slid));
Dict_mask = sum(DictsA,1);
DictsA(:,Dict_mask<=0.4)=[]; 
vecOfMeans = mean(DictsA);
DictsA = DictsA-ones(size(DictsA,1),1)*vecOfMeans;
dictA_Atom  = double(DictsA(:,1:50:end)); clear DictsA
DictsB= single(im2colstep(double(Noise), par.b, par.slid));
DictsB(:,Dict_mask<=0.4)=[];
vecOfMeans = mean(DictsB);
DictsB = DictsB-ones(size(DictsB,1),1)*vecOfMeans;
dictB_Atom  = double(DictsB(:,1:50:end));clear DictsB
ratioA=spdiag(double(1./sqrt(sum(dictA_Atom.*dictA_Atom))));  
dictA_Atom = dictA_Atom*ratioA;%%
ratioB=spdiag(double(1./sqrt(sum(dictB_Atom.*dictB_Atom))));
dictB_Atom = dictB_Atom*ratioB;
dictA_Atom =dictA_Atom(:,1:4:end);
dictB_Atom =dictB_Atom(:,1:4:end);
% save('C:\Users\Liuhanxi\Desktop\Recon_matlab\MAT\DFR_64_4_SIRT','dictA','dictB')

par.b = par.b(1);
dictimg = showdict(dictA_Atom(1:par.b*par.b,:),[1 1]*par.b,round(sqrt(size(dictA_Atom,2))),round(sqrt(size(dictA_Atom,2))),'lines','highcontrast');
figure, imshow(dictimg);
dictimg = showdict(dictB_Atom(1:par.b*par.b,:),[1 1]*par.b,round(sqrt(size(dictB_Atom,2))),round(sqrt(size(dictB_Atom,2))),'lines','highcontrast');
figure, imshow(dictimg);


