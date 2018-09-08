function mask = spm_load_mask(V)
% Loads the mask for registration
% Adapted from:
% ____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_load_priors8.m 3168 2009-05-29 20:52:52Z john $

tiny = 1e-5;

deg = 1;
if ~isstruct(V), V  = spm_vol(V); end;
% spm_check_orientations(V);

mask.V = V;
mask.M = mask.V(1).mat;

Kb = numel(mask.V);
mask.dat = cell(Kb,1);
for k1=1:(Kb),
    mask.dat{k1} = zeros(mask.V(1).dim(1:3));
end;

spm_progress_bar('Init',mask.V(1).dim(3),'Loading priors','Planes loaded');
for i=1:mask.V(1).dim(3)
    M         = spm_matrix([0 0 i]);
%     s         = zeros(mask.V(1).dim(1:2));
    for k1=1:Kb
        tmp                = spm_slice_vol(mask.V(k1),M,mask.V(1).dim(1:2),0);
        mask.dat{k1}(:,:,i) = max(min(tmp,1),0);
%         s                  = s + tmp;
    end;
%     t = s>1;
%     if any(t)
%         for k1=1:Kb
%             tmp           = mask.dat{k1}(:,:,i);
%             tmp(t)        = tmp(t)./s(t); % normalize the loaded mask
%             % emask.nii will be normalized here. Others (bmask, cmask, cmaskthresh) were normalized before loading.
%             mask.dat{k1}(:,:,i) = tmp;
%         end;
%     end;
    spm_progress_bar('Set',i);
end;
% mask.bg1 = zeros(Kb,1);
for k1=1:Kb,
%     mask.bg1(k1)  = mean(mean(mask.dat{k1}(:,:,1)));
%     mask.bg2(k1)  = mean(mean(mask.dat{k1}(:,:,end)));
    mask.dat{k1} = spm_bsplinc(log(mask.dat{k1}+tiny),[deg deg deg  0 0 0]); % get B-spline coefficients for subsequent sampling using B-spline interpolation
end;
% mask.tiny = tiny;
mask.deg  = deg+1;

spm_progress_bar('Clear');
return;

