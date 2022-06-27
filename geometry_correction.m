% function proj_geom = geometry_correction(pMatrix_file, nview,nu, proj_geom)
function proj_geom = geometry_correction(pMatrix_file, nu, proj_geom)
nview = numel(proj_geom.ProjectionAngles);
du = proj_geom.DetectorSpacingX;
dv = proj_geom.DetectorSpacingY;
nv = proj_geom.DetectorRowCount;
nu = nu;

fid = fopen(pMatrix_file,'rb');
pMatrix = fread(fid, [12,nview], 'float');
pMatrix = reshape(pMatrix, [3,4,nview]);
pOffsetU = -(squeeze((pMatrix(1, 4, :) ./ pMatrix(3, 4, :) - (nu/2.0+0.5))));
pOffsetV = -(squeeze((pMatrix(2, 4, :) ./ pMatrix(3, 4, :) - (nv/2.0+0.5))));
fclose(fid);

proj_geom =astra_geom_2vec(proj_geom);
for i=1:nview
    proj_geom.Vectors(i, 4:6) = proj_geom.Vectors(i, 4:6)...
        + pOffsetU(i) * proj_geom.Vectors(i, 7:9)...
        + pOffsetV(i) * proj_geom.Vectors(i, 10:12);
end

end

