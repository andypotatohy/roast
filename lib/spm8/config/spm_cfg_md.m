function md = spm_cfg_md
% SPM Configuration file for 'mkdir'
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_md.m 4907 2012-09-06 19:33:21Z guillaume $

% ----------------------------------------------------------------------
% basedir Select a base directory
% ----------------------------------------------------------------------
basedir         = cfg_files;
basedir.tag     = 'basedir';
basedir.name    = 'Select a base directory';
basedir.help    = {'Select a base directory.'};
basedir.filter  = 'dir';
basedir.ufilter = '.*';
basedir.num     = [1 1];

% ----------------------------------------------------------------------
% name Enter a directory name
% ----------------------------------------------------------------------
name            = cfg_entry;
name.tag        = 'name';
name.name       = 'Enter a directory name';
name.help       = {'Enter a directory name'};
name.strtype    = 's';
name.num        = [1 Inf];

% ----------------------------------------------------------------------
% md Make Directory (Deprecated)
% ----------------------------------------------------------------------
md              = cfg_exbranch;
md.tag          = 'md';
md.name         = 'Make Directory (DEPRECATED)';
md.val          = {basedir name};
md.help         = {
    'This module is DEPRECATED and has been moved to BasicIO.'
    'Jobs which are ready to run may continue using it, but the module inputs can not be changed via GUI.'
    'Please switch to the BasicIO module instead.'
    'This module will be REMOVED in the next major release of SPM.'
    ''
    'This facility allows programming a directory change.'
}';
md.prog         = @my_mkdir;
md.hidden       = true;

%=======================================================================
function my_mkdir(varargin)
warning('"spm.util.md" is DEPRECATED and will be REMOVED in the next major release of SPM. Use BasicIO instead.');
job = varargin{1};
if ~isempty(job.basedir) && ~isempty(job.name)
    mkdir(job.basedir{:},job.name);
end
