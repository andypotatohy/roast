function ef = getEF(fm,coords,mon,mask)

Ncoords = size(coords,1);
[xi,yi,zi] = ndgrid(1:size(mask,1),1:size(mask,2),1:size(mask,3));

ef = zeros([size(mask) 3]);

xopt = fm*mon;

Fx = TriScatteredInterp(coords, xopt(1:Ncoords));
ef(:,:,:,1) = Fx(xi,yi,zi);
Fy = TriScatteredInterp(coords, xopt(Ncoords+1:2*Ncoords));
ef(:,:,:,2) = Fy(xi,yi,zi);
Fz = TriScatteredInterp(coords, xopt(2*Ncoords+1:3*Ncoords));
ef(:,:,:,3) = Fz(xi,yi,zi);
