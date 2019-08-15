function [efVol,ef] = getEF(fm,coords,mon,dim)

Ncoords = size(coords,1);
[xi,yi,zi] = ndgrid(1:dim(1),1:dim(2),1:dim(3));

efVol = zeros([dim 3]);

ef = fm*mon;

Fx = TriScatteredInterp(coords, ef(1:Ncoords));
efVol(:,:,:,1) = Fx(xi,yi,zi);
Fy = TriScatteredInterp(coords, ef(Ncoords+1:2*Ncoords));
efVol(:,:,:,2) = Fy(xi,yi,zi);
Fz = TriScatteredInterp(coords, ef(2*Ncoords+1:3*Ncoords));
efVol(:,:,:,3) = Fz(xi,yi,zi);
