
function  Rec = SART_Cone_dict(xint,proj_geom,vol_geom,proj,Norimg,IM,IW,lambda,par)
Er          = 1e-10;
Rec         = xint;
img_diff    = zeros(size(xint),'single'); 
id          = Norimg > Er;
for i = 1: par.niter
    [proj_id, ARec]    =  astra_create_sino3d_cuda(Rec, proj_geom, vol_geom);  
    astra_mex_data3d('delete', proj_id);
    proj_diff          =  ARec-proj;
    [vol_id,tmp]       =  astra_create_backprojection3d_cuda(proj_diff, proj_geom, vol_geom);  
    astra_mex_data3d('delete', vol_id);
             tmp2           = lambda*(IW.*Rec-IM);
             img_diff(id)   = (tmp(id)+tmp2(id))./(Norimg(id)+lambda*IW(id));
             Rec            = Rec-img_diff.*par.mask;
             Rec            = max(Rec,0);
%            Rec            = min(Rec,par.pixmax);
end
end