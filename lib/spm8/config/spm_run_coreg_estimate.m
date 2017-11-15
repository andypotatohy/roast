function out = spm_run_coreg_estimate(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg_estimate.m 4380 2011-07-05 11:27:12Z volkmar $

job = varargin{1};
if isempty(job.other{1})
    job.other = {};
end

x  = spm_coreg(char(job.ref), char(job.source), job.eoptions);

M  = spm_matrix(x);
PO = [job.source(:); job.other(:)];
MM = zeros(4,4,numel(PO));
for j=1:numel(PO)
    MM(:,:,j) = spm_get_space(PO{j});
end
for j=1:numel(PO)
    spm_get_space(PO{j}, M\MM(:,:,j));
end

% % HACK ONLY
% job.source{1} = '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,1';
% job.other{1} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,2';
% job.other{2} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,3';
% job.other{3} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,4';
% job.other{4} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,5';
% job.other{5} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,6';
% 
% PO = [job.source(:); job.other(:)];
% MM = repmat(MM,[1,1,6]);
% 
% for j=1:numel(PO)
%     spm_get_space(PO{j}, M\MM(:,:,j));
% end
% % HACK ONLY
% % ANDY 2016-08-24

out.cfiles = PO;
out.M      = M;