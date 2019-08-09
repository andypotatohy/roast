function roast_target(subj,simTag,targetCoord,varargin)

fprintf('\n\n');
disp('=============================================================')
disp('ROAST is an aggregated work by Yu (Andy) Huang licensed under')
disp('General Public License version 3 or later. It''s supported by')
disp('both NIH grants and Soterix Medical Inc.')
disp('=============================================================')

addpath(genpath([fileparts(which(mfilename)) filesep 'lib/']));

fprintf('\n\n');
disp('======================================================')
disp('CHECKING INPUTS...')
disp('======================================================')
fprintf('\n');

% warning('on');

% check subject name
if nargin<1 || isempty(subj)
    subj = 'example/MNI152_T1_1mm.nii';
end

if strcmpi(subj,'nyhead')
    subj = 'example/nyhead.nii';
end

% if ~strcmpi(subj,'example/nyhead.nii') && ~exist(subj,'file')
%     error(['The subject MRI you provided ' subj ' does not exist.']);
% end

% check simulation tag
if nargin<2 || isempty(simTag)
    error(['Please provide a valid simulation tag for Subject ' subj ', so that roast_target can locate the corresponding lead field.']);
end

[dirname,baseFilename,ext] = fileparts(subj);
if ~exist([dirname filesep baseFilename '_' simTag '_options.mat'],'file')
    error(['Option file not found. Simulation ' simTag ' may never be run. Please run it first.']);
else
    load([dirname filesep baseFilename '_' simTag '_options.mat'],'opt');
end

if ~strcmp(opt.configTxt,'leadFieldGeneration')
    error(['Simulation ' simTag ' was NOT run for generating the lead field. roast_target() cannot work without a proper lead field. Please run roast with ''leadField'' as the recipe first.']);
end

% check target coordinates
if nargin<3 || isempty(targetCoord)
    targetCoord = [-48 -8 50]; % defaults to left primary motor cortex
end

if size(targetCoord,2) ~= 3
    error('Unrecognized format of target coordinates. Please enter as [x y z].');
end

% take in user-specified options
if mod(length(varargin),2)~=0
    error('Unrecognized format of options. Please enter as property-value pair.');
end

indArg = 1;
while indArg <= length(varargin)
    switch lower(varargin{indArg})
        case 'coordtype'
            coordType = varargin{indArg+1};
            indArg = indArg+2;
        case 'opttype'
            optType = varargin{indArg+1};
            indArg = indArg+2;
        case 'orient'
            orient = varargin{indArg+1};
            indArg = indArg+2;
        case 'elecnum'
            elecNum = varargin{indArg+1};
            indArg = indArg+2;
        case 'targetradius'
            targetRadius = varargin{indArg+1};
            indArg = indArg+2;
        case 'k'
            k = varargin{indArg+1};
            indArg = indArg+2;
        case 'a'
            a = varargin{indArg+1};
            indArg = indArg+2;        
        otherwise
            error('Supported options are: ''coordType'', ''optType'', ''orient'', ''elecNum'', ''targetRadius'', ''k'' and ''a''.');
    end
end

% set up defaults and check on option conflicts
if ~exist('coordType','var')
    coordType = 'mni';
else
    if ~any(strcmpi(coordType,{'mni','voxel'}))
        error('Please enter either ''mni'' or ''voxel'' for option ''coordType''.');
    end
end

if ~exist('optType','var')
    optType = 'max-l1';
else
    if ~any(strcmpi(optType,{'unconstrained-wls','wls-l1','wls-l1per','wls-l1penalty','unconstrained-lcmv','lcmv-l1','lcmv-l1per','max-l1','max-l1per'}))
        error('Supported targeting optimizations are: ''unconstrained-wls'',''wls-l1'',''wls-l1per'',''wls-l1penalty'',''unconstrained-lcmv'',''lcmv-l1'',''lcmv-l1per'',''max-l1'' and ''max-l1per''.');
    end
end

if ~exist('orient','var')
    orient = 'radial-in';
else
    if ~any(strcmpi(orient,{'radial-in','radial-out','right','left','anterior','posterior','right-anterior','right-posterior','left-anterior','left-posterior','optimal'}))
        error('Supported orientations of electric field at the target are: ''radial-in'',''radial-out'',''right'',''left'',''anterior'',''posterior'',''right-anterior'',''right-posterior'',''left-anterior'',''left-posterior'',''optimal''.');
    end
end

if ~exist('elecNum','var')
    if strcmpi(optType,'max-l1per')
        elecNum = 4;
    else
        elecNum = [];
    end
else
    if strcmpi(optType,'max-l1per')
        if elecNum<4 || mod(elecNum,1)~=0
            error('Unrecognized option value. Please enter positive integer value of at least 4 for option ''elecNum''.');
        end
    else
        warning('Option ''elecNum'' will be ignored for optimization type other than ''max-l1per''.');
        elecNum = [];
    end
end

if ~exist('targetRadius','var')
    targetRadius = 3;
else
    if targetRadius<=0 || mod(targetRadius,1)~=0
        error('Unrecognized option value. Please enter positive integer value for option ''targetRadius''.');
    end
end

if ~exist('k','var')
    k = 0.8;
else
    if k<=0 || k>=1
        error('Unrecognized option value. Please enter the ''k'' value between 0 and 1.');
    end
    warning('You''re changing the advanced options of ROAST-TARGET. Unless you know what you''re doing, please keep the ''k'' value default.');
end

if ~exist('a','var')
    a = 0.5;
else
    if a<=0 || a>=1
        error('Unrecognized option value. Please enter the ''a'' value between 0 and 1.');
    end
    warning('You''re changing the advanced options of ROAST-TARGET. Unless you know what you''re doing, please keep the ''a'' value default.');
end