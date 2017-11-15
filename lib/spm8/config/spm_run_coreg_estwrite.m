function out = spm_run_coreg_estwrite(varargin)
% SPM job execution function
% takes a harvested job data structure and call SPM functions to perform
% computations on the data.
% Input:
% job    - harvested job data structure (see matlabbatch help)
% Output:
% out    - computation results, usually a struct variable.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_coreg_estwrite.m 4380 2011-07-05 11:27:12Z volkmar $

job = varargin{1};
if isempty(job.other{1})
    job.other = {};
end

x  = spm_coreg(char(job.ref), char(job.source),job.eoptions);

M  = spm_matrix(x);
PO = [job.source(:); job.other(:)];
MM = zeros(4,4,numel(PO));
for j=1:numel(PO),
    MM(:,:,j) = spm_get_space(PO{j});
end

% % hack ONLY for project of NYUpatient
% P1 = job.ref(:);
% sv2rv = inv(spm_get_space(P1{1}))*inv(M)*MM;
% % sv2rw =                             inv(M)*MM;
% % source voxel to ref world space % ANDY 2016-03-18
% % but this actually did NOT work, because ScanIP generates model in
% % pseudo-world space, ie., only considers scaling, no shifting or shearing.
% % So cannot use formal mapping to world space to locate the electrodes, and
% % have to use hack (just hack the scaling from voxel space coordinates).
% % ANDY 2016-03-18
% save('sv2rv.mat','sv2rv');
% % save('sv2rw.mat','sv2rw');
% % ANDY 2015-06-02
% 
% % % HACK ONLY
% % job.source{1} = '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,1';
% % job.other{1} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,2';
% % job.other{2} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,3';
% % job.other{3} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,4';
% % job.other{4} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,5';
% % job.other{5} =  '/home/andy/projects/modelValidation/data/patientFromNYU/eTPMn.nii,6';
% % 
% % PO = [job.source(:); job.other(:)];
% % MM = repmat(MM,[1,1,6]);
% % % HACK ONLY
% % % ANDY 2016-08-24

for j=1:numel(PO),
    spm_get_space(PO{j}, M\MM(:,:,j));
end

P            = char(job.ref{:},job.source{:},job.other{:});
flags.mask   = job.roptions.mask;
flags.mean   = 0;
flags.interp = job.roptions.interp;
flags.which  = 1;
flags.wrap   = job.roptions.wrap;
flags.prefix = job.roptions.prefix;

spm_reslice(P,flags);
% spm_reslice_hack(P,flags); % ANDY 2014-02-14

out.cfiles = PO;
out.M      = M;
out.rfiles = cell(size(out.cfiles));
for i=1:numel(out.cfiles),
    [pth,nam,ext,num] = spm_fileparts(out.cfiles{i});
    out.rfiles{i} = fullfile(pth,[job.roptions.prefix, nam, ext, num]);
end;
return;

