function cdir = spm_cfg_cdir
% SPM Configuration file for 'cd'
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_cdir.m 4907 2012-09-06 19:33:21Z guillaume $

rev = '$Rev: 4907 $';
% ---------------------------------------------------------------------
% directory Select a directory
% ---------------------------------------------------------------------
directory         = cfg_files;
directory.tag     = 'directory';
directory.name    = 'Select a directory';
directory.help    = {'Select a directory to change to.'};
directory.filter = 'dir';
directory.ufilter = '.*';
directory.num     = [1 1];
% ---------------------------------------------------------------------
% cdir Change Directory (Deprecated)
% ---------------------------------------------------------------------
cdir         = cfg_exbranch;
cdir.tag     = 'cdir';
cdir.name    = 'Change Directory (DEPRECATED)';
cdir.val     = {directory };
cdir.help    = {
    'This module is DEPRECATED and has been moved to BasicIO.'
    'Jobs which are ready to run may continue using it, but the module inputs can not be changed via GUI.'
    'Please switch to the BasicIO module instead.'
    'This module will be REMOVED in the next major release of SPM.'
    ''
    'This facility allows programming a directory change.'
}';
cdir.prog = @my_job_cd;
cdir.hidden = true;
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function my_job_cd(varargin)
% job can be a job structure or the directory to change to.
warning('"spm.util.cdir" is DEPRECATED and will be REMOVED in the next major release of SPM. Use BasicIO instead.');
job = varargin{1};
if isstruct(job)
    jobDir = job.directory;
else
    jobDir = job;
end
if ~isempty(jobDir)
    cd(char(jobDir));
    fprintf('New working directory: %s\n', char(jobDir));
end
