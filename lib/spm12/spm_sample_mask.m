function s = spm_sample_mask(mask,x1,x2,x3)
% Sample mask to reference space
% Adapted from:
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_sample_priors8.m 3168 2009-05-29 20:52:52Z john $

deg  = mask.deg;
% tiny = mask.tiny;

d  = size(mask.dat{1});
dx = size(x1);
Kb = numel(mask.dat);
s  = cell(1,Kb);
msk1 = x1>=1 & x1<=d(1) & x2>=1 & x2<=d(2) & x3>=1 & x3<=d(3);
% msk2 = x3<1;
x1 = x1(msk1);
x2 = x2(msk1);
x3 = x3(msk1);
% if nargout<=1,
%     tot = zeros(dx);
    for k=1:Kb,
        a    = spm_bsplins(mask.dat{k},x1,x2,x3,[deg deg deg  0 0 0]); % sampling mask using B-spline interpolation
                                                % deg is NOT same as used by spm_bsplinc???
%         s{k} = ones(dx)*mask.bg2(k);
s{k} = ones(dx)*0;
        s{k}(msk1) = exp(a);
%         s{k}(msk2) = mask.bg1(k);
%         tot  = tot + s{k};
    end
%     msk      = ~isfinite(tot);
%     tot(msk) = 1;
%     for k=1:Kb,
%         s{k}(msk) = mask.bg2(k);
%         s{k}      = s{k}./tot; % normalize the loaded mask
%     end
% else
%     ds1 = cell(1,Kb);
%     ds2 = cell(1,Kb);
%     ds3 = cell(1,Kb);
%     tot = zeros(dx);
%     for k=1:Kb,
%         [a,da1,da2,da3] = spm_bsplins(mask.dat{k},x1,x2,x3,[deg deg deg  0 0 0]);
%         if k==Kb, s{k} = ones(dx); else s{k} = zeros(dx)+tiny; end
%         s{k} = ones(dx)*mask.bg2(k);
%         s{k}(msk1) = exp(a);
%         s{k}(msk2) = mask.bg1(k);
%         tot    = tot + s{k};
%         ds1{k} = zeros(dx); ds1{k}(msk1) = da1;
%         ds2{k} = zeros(dx); ds2{k}(msk1) = da2;
%         ds3{k} = zeros(dx); ds3{k}(msk1) = da3;
%     end
%     msk      = find(~isfinite(tot));
%     tot(msk) = 1;
%     da1      = zeros(dx);
%     da2      = zeros(dx);
%     da3      = zeros(dx);
%     for k=1:Kb,
%          s{k}(msk) = mask.bg1(k);
%          s{k}      = s{k}./tot;
%          da1       = da1 + s{k}.*ds1{k};
%          da2       = da2 + s{k}.*ds2{k};
%          da3       = da3 + s{k}.*ds3{k};
%     end
%     for k=1:Kb,
%         ds1{k} = s{k}.*(ds1{k} - da1); ds1{k}(msk) = 0;
%         ds2{k} = s{k}.*(ds2{k} - da2); ds2{k}(msk) = 0;
%         ds3{k} = s{k}.*(ds3{k} - da3); ds3{k}(msk) = 0;
%     end
% end;

