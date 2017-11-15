function [xy,label] = spm_eeg_layout3D(sens, modality)
% Wrapper function to a fieldtrip function to project 3D locations 
% onto a 2D plane. 
% FORMAT [xy,label] = spm_eeg_project3D(sens, modality)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_layout3D.m 5141 2012-12-20 15:54:55Z vladimir $


switch modality
    case 'EEG'
        xy = sens.chanpos;
        label = sens.label;
    case 'MEG'
        cfg = [];
        cfg.style = '3d';
        cfg.rotate = 0;
        cfg.grad = sens;
        cfg.showcallinfo = 'no';

        lay = ft_prepare_layout(cfg);

        [sel1, sel2] = spm_match_str(sens.label, lay.label);

        label = lay.label(:);
        label = label(sel2, 1);
        xy = lay.pos(sel2, :);
end

