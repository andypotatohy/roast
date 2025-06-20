function roast(subj,recipe,varargin)
% roast(subj,recipe,varargin)
%
% Main function of ROAST.
% 
% Please refer to the README.md on the github repo for better formated
% documentations: https://github.com/andypotatohy/roast
% 
% If you use ROAST in your research, please cite these:
% 
% Huang, Y., Datta, A., Bikson, M., Parra, L.C., Realistic vOlumetric-Approach
% to Simulate Transcranial Electric Stimulation -- ROAST -- a fully automated
% open-source pipeline, Journal of Neural Engineering, Vol. 16, No. 5, 2019 (prefered reference)
% 
% Huang, Y., Datta, A., Bikson, M., Parra, L.C., ROAST: an open-source,
% fully-automated, Realistic vOlumetric-Approach-based Simulator for TES,
% Proceedings of the 40th Annual International Conference of the IEEE Engineering
% in Medicine and Biology Society, Honolulu, HI, July 2018
% 
% If you use New York head to run simulation, please also cite the following:
% Huang, Y., Parra, L.C., Haufe, S.,2016. The New York Head - A precise
% standardized volume conductor model for EEG source localization and tES
% targeting, NeuroImage,140, 150-162
% 
% If you also use the targeting feature (`roast_target`), please cite these:
% 
% Dmochowski, J.P., Datta, A., Bikson, M., Su, Y., Parra, L.C., Optimized 
% multi-electrode stimulation increases focality and intensity at target,
% Journal of Neural Engineering 8 (4), 046011, 2011
% 
% Dmochowski, J.P., Datta, A., Huang, Y., Richardson, J.D., Bikson, M.,
% Fridriksson, J., Parra, L.C., Targeted transcranial direct current stimulation 
% for rehabilitation after stroke, NeuroImage, 75, 12-19, 2013
% 
% Huang, Y., Thomas, C., Datta, A., Parra, L.C., Optimized tDCS for Targeting
% Multiple Brain Regions: An Integrated Implementation. Proceedings of the 40th
% Annual International Conference of the IEEE Engineering in Medicine and Biology
% Society, Honolulu, HI, July 2018, 3545-3548
% 
% If you also use the Multiaxial for segmentation by turning on the `multiaxial`
% option, please cite this:
%
% Full-Head Segmentation of MRI with Abnormal Brain Anatomy: Model and Data Release
% Andrew M Birnbaum, Adam Buchwald, Peter Turkeltaub, Adam Jacks, Yu Huang, Abhisheck Datta, Lucas C Parra, Lukas A Hirsch
% https://arxiv.org/abs/2501.18716 (Jan 2025).
% 
% ROAST was supported by the NIH through grants R01MH111896, R01MH111439, 
% R01NS095123, R44NS092144, R41NS076123, and by Soterix Medical Inc.
% 
% General Public License version 3 or later. See LICENSE.md for details.
% 
% This software uses free packages from the Internet, except Matlab, which
% is a proprietary software by the MathWorks. You need a valid Matlab license
% to run this software.
% 
% ROAST is considered as an "aggregate" rather than "derived work", based on
% the definitions in GPL FAQ. The ROAST license only applies to the scripts,
% documentation and the individual MRI data under example/ folder in this 
% package and excludes those programs stored in the lib/ directory. The software 
% under lib/ follow their respective licenses. This software is only intended
% for non-commercial use.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% September 2019
%
% (c) Multiaxial segmentation developed by Lukas Hirsch, and integrated into
% ROAST by Andrew Birnbaum
% April 2024

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

if ~strcmpi(subj,'example/nyhead.nii') && ~exist(subj,'file')
    error(['The subject MRI you provided ' subj ' does not exist.']);
end

if nargin<2 || isempty(recipe)
    recipe = {'Fp1',1,'P4',-1};
end

% take in user-specified options
if mod(length(varargin),2)~=0
    error('Unrecognized format of options. Please enter as property-value pair.');
end

indArg = 1;
while indArg <= length(varargin)
    switch lower(varargin{indArg})
        case 'captype'
            capType = varargin{indArg+1};
            indArg = indArg+2;
        case 'electype'
            elecType = varargin{indArg+1};
            indArg = indArg+2;
        case 'elecsize'
            elecSize = varargin{indArg+1};
            indArg = indArg+2;
        case 'elecori'
            elecOri = varargin{indArg+1};
            indArg = indArg+2;
        case 't2'
            T2 = varargin{indArg+1};
            indArg = indArg+2;
        case 'multiaxial'
            multiaxial = varargin{indArg+1};
            indArg = indArg+2;
        case 'manual_gui'
            manual_gui = varargin{indArg+1};
            indArg = indArg+2;
        case 'meshoptions'
            meshOpt = varargin{indArg+1};
            indArg = indArg+2;
        case 'conductivities'
            conductivities = varargin{indArg+1};
            indArg = indArg+2;
        case 'simulationtag'
            simTag = varargin{indArg+1};
            indArg = indArg+2;
        case 'resampling'
            doResamp = varargin{indArg+1};
            indArg = indArg+2;
        case 'zeropadding'
            paddingAmt = varargin{indArg+1};
            indArg = indArg+2;
        otherwise
            error('Supported options are: ''capType'', ''elecType'', ''elecSize'', ''elecOri'', ''T2'', ''multiaxial'',''manual_gui'',''meshOptions'',''conductivities'', ''simulationTag'', ''resampling'', and ''zeroPadding''.');
    end
end

if any(~strcmpi(recipe,'leadfield'))
    
    % check recipe syntax
    if mod(length(recipe),2)~=0
        error('Unrecognized format of your recipe. Please enter as electrodeName-injectedCurrent pair.');
    end
    
    elecName = (recipe(1:2:end-1))';
    injectCurrent = (cell2mat(recipe(2:2:end)))';
    if abs(sum(injectCurrent))>eps
        error('Electric currents going in and out of the head not balanced. Please make sure they sum to 0.');
    end
    
    % set up defaults and check on option conflicts
    if ~exist('capType','var')
        capType = '1010';
    else
        if ~any(strcmpi(capType,{'1020','1010','1005','biosemi','egi'}))
            error('Supported cap types are: ''1020'', ''1010'', ''1005'', ''BioSemi'' and ''EGI''.');
        end
    end
    
    if ~exist('elecType','var')
        elecType = 'disc';
    else
        if ~iscellstr(elecType)
            if ~any(strcmpi(elecType,{'disc','pad','ring'}))
                error('Supported electrodes are: ''disc'', ''pad'' and ''ring''.');
            end
        else
            if length(elecType)~=length(elecName)
                error('You want to place more than 1 type of electrodes, but did not tell ROAST which type for each electrode. Please provide the type for each electrode respectively, as the value for option ''elecType'', in a cell array of length equals to the number of electrodes to be placed.');
            end
            for i=1:length(elecType)
                if ~any(strcmpi(elecType{i},{'disc','pad','ring'}))
                    error('Supported electrodes are: ''disc'', ''pad'' and ''ring''.');
                end
            end
        end
    end
    
    if ~exist('elecSize','var')
        if ~iscellstr(elecType)
            switch lower(elecType)
                case {'disc'}
                    elecSize = [6 2];
                case {'pad'}
                    elecSize = [50 30 3];
                case {'ring'}
                    elecSize = [4 6 2];
            end
        else
            elecSize = cell(1,length(elecType));
            for i=1:length(elecSize)
                switch lower(elecType{i})
                    case {'disc'}
                        elecSize{i} = [6 2];
                    case {'pad'}
                        elecSize{i} = [50 30 3];
                    case {'ring'}
                        elecSize{i} = [4 6 2];
                end
            end
        end
    else
        if ~iscellstr(elecType)
            if iscell(elecSize)
                warning('Looks like you''re placing only 1 type of electrodes. ROAST will only use the 1st entry of the cell array of ''elecSize''. If this is not what you want and you meant differect sizes for different electrodes of the same type, just enter ''elecSize'' option as an N-by-2 or N-by-3 matrix, where N is number of electrodes to be placed.');
                elecSize = elecSize{1};
            end
            if any(elecSize(:)<=0)
                error('Please enter non-negative values for electrode size.');
            end
            if size(elecSize,2)~=2 && size(elecSize,2)~=3
                error('Unrecognized electrode sizes. Please specify as [radius height] for disc, [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
            end
            if size(elecSize,1)>1 && size(elecSize,1)~=length(elecName)
                error('You want different sizes for each electrode. Please tell ROAST the size for each electrode respectively, in a N-row matrix, where N is the number of electrodes to be placed.');
            end
            if strcmpi(elecType,'disc') && size(elecSize,2)==3
                error('Redundant size info for Disc electrodes. Please enter as [radius height]');
                %             elecSize = elecSize(:,1:2);
            end
            if any(strcmpi(elecType,{'pad','ring'})) && size(elecSize,2)==2
                error('Insufficient size info for Pad or Ring electrodes. Please specify as [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
            end
            if strcmpi(elecType,'pad') && any(elecSize(:,1) < elecSize(:,2))
                error('For Pad electrodes, the width of the pad should not be bigger than its length. Please enter as [length width height]');
            end
            if strcmpi(elecType,'pad') && any(elecSize(:,3) < 3)
                error('For Pad electrodes, the thickness should at least be 3 mm.');
            end
            if strcmpi(elecType,'pad') && any(elecSize(:) > 80)
                warning('You''re placing large pad electrodes (one of its dimensions is bigger than 8 cm). For large pads, the size will not be exact in the model because they will be bent to fit the scalp surface.');
            end
            if strcmpi(elecType,'ring') && any(elecSize(:,1) >= elecSize(:,2))
                error('For Ring electrodes, the inner radius should be smaller than outter radius. Please enter as [innerRadius outterRadius height]');
            end
        else
            if ~iscell(elecSize)
                error('You want to place at least 2 types of electrodes, but only provided size info for 1 type. Please provide complete size info for all types of electrodes in a cell array as the value for option ''elecSize'', or just use defaults by not specifying ''elecSize'' option.');
            end
            if length(elecSize)~=length(elecType)
                error('You want to place more than 1 type of electrodes. Please tell ROAST the size for each electrode respectively, as the value for option ''elecSize'', in a cell array of length equals to the number of electrodes to be placed.');
            end
            for i=1:length(elecSize)
                if isempty(elecSize{i})
                    switch lower(elecType{i})
                        case {'disc'}
                            elecSize{i} = [6 2];
                        case {'pad'}
                            elecSize{i} = [50 30 3];
                        case {'ring'}
                            elecSize{i} = [4 6 2];
                    end
                else
                    if any(elecSize{i}(:)<=0)
                        error('Please enter non-negative values for electrode size.');
                    end
                    if size(elecSize{i},2)~=2 && size(elecSize{i},2)~=3
                        error('Unrecognized electrode sizes. Please specify as [radius height] for disc, [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
                    end
                    if size(elecSize{i},1)>1
                        error('You''re placing more than 1 type of electrodes. Please put size info for each electrode as a 1-row vector in a cell array for option ''elecSize''.');
                    end
                    if strcmpi(elecType{i},'disc') && size(elecSize{i},2)==3
                        error('Redundant size info for Disc electrodes. Please enter as [radius height]');
                        %                     elecSize{i} = elecSize{i}(:,1:2);
                    end
                    if any(strcmpi(elecType{i},{'pad','ring'})) && size(elecSize{i},2)==2
                        error('Insufficient size info for Pad or Ring electrodes. Please specify as [length width height] for pad, and [innerRadius outterRadius height] for ring electrode.');
                    end
                    if strcmpi(elecType{i},'pad') && any(elecSize{i}(:,1) < elecSize{i}(:,2))
                        error('For Pad electrodes, the width of the pad should not be bigger than its length. Please enter as [length width height]');
                    end
                    if strcmpi(elecType{i},'pad') && any(elecSize{i}(:,3) < 3)
                        error('For Pad electrodes, the thickness should at least be 3 mm.');
                    end
                    if strcmpi(elecType{i},'pad') && any(elecSize{i}(:) > 80)
                        warning('You''re placing large pad electrodes (one of its dimensions is bigger than 8 cm). For large pads, the size will not be exact in the model because they will be bent to fit the scalp surface.');
                    end
                    if strcmpi(elecType{i},'ring') && any(elecSize{i}(:,1) >= elecSize{i}(:,2))
                        error('For Ring electrodes, the inner radius should be smaller than outter radius. Please enter as [innerRadius outterRadius height]');
                    end
                end
            end
        end
    end
    
    if ~exist('elecOri','var')
        if ~iscellstr(elecType)
            if strcmpi(elecType,'pad')
                elecOri = 'lr';
            else
                elecOri = [];
            end
        else
            elecOri = cell(1,length(elecType));
            for i=1:length(elecOri)
                if strcmpi(elecType{i},'pad')
                    elecOri{i} = 'lr';
                else
                    elecOri{i} = [];
                end
            end
        end
    else
        if ~iscellstr(elecType)
            if ~strcmpi(elecType,'pad')
                warning('You''re not placing pad electrodes; customized orientation options will be ignored.');
                elecOri = [];
            else
                if iscell(elecOri)
                    allChar = 1;
                    for i=1:length(elecOri)
                        if ~ischar(elecOri{i})
                            warning('Looks like you''re only placing pad electrodes. ROAST will only use the 1st entry of the cell array of ''elecOri''. If this is not what you want and you meant differect orientations for different pad electrodes, just enter ''elecOri'' option as an N-by-3 matrix, or as a cell array of length N (put ''lr'', ''ap'', or ''si'' into the cell element), where N is number of pad electrodes to be placed.');
                            elecOri = elecOri{1};
                            allChar = 0;
                            break;
                        end
                    end
                    if allChar && length(elecOri)~=length(elecName)
                        error('You want different orientations for each pad electrode by using pre-defined keywords in a cell array. Please make sure the cell array has a length equal to the number of pad electrodes.');
                    end
                end
                if ~iscell(elecOri)
                    if ischar(elecOri)
                        if ~any(strcmpi(elecOri,{'lr','ap','si'}))
                            error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                        end
                    else
                        if size(elecOri,2)~=3
                            error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                        end
                        if size(elecOri,1)>1 && size(elecOri,1)~=length(elecName)
                            error('You want different orientations for each pad electrode. Please tell ROAST the orientation for each pad respectively, in a N-by-3 matrix, where N is the number of pads to be placed.');
                        end
                    end
                end
            end
        else
            if ~iscell(elecOri)
                elecOri0 = elecOri;
                elecOri = cell(1,length(elecType));
                if ischar(elecOri0)
                    if ~any(strcmpi(elecOri0,{'lr','ap','si'}))
                        error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                    end
                    for i=1:length(elecType)
                        if strcmpi(elecType{i},'pad')
                            elecOri{i} = elecOri0;
                        else
                            elecOri{i} = [];
                        end
                    end
                else
                    if size(elecOri0,2)~=3
                        error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                    end
                    numPad = 0;
                    for i=1:length(elecType)
                        if strcmpi(elecType{i},'pad')
                            numPad = numPad+1;
                        end
                    end
                    if size(elecOri0,1)>1
                        if size(elecOri0,1)~=numPad
                            error('You want different orientations for each pad electrode. Please tell ROAST the orientation for each pad respectively, in a N-by-3 matrix, where N is the number of pads to be placed.');
                        end
                    else
                        elecOri0 = repmat(elecOri0,numPad,1);
                    end
                    i0=1;
                    for i=1:length(elecType)
                        if strcmpi(elecType{i},'pad')
                            elecOri{i} = elecOri0(i0,:);
                            i0 = i0+1;
                        else
                            elecOri{i} = [];
                        end
                    end
                end
            else
                if length(elecOri)~=length(elecType)
                    error('You want to place another type of electrodes aside from pad. Please tell ROAST the orienation for each electrode respectively, as the value for option ''elecOri'', in a cell array of length equals to the number of electrodes to be placed (put [] for non-pad electrodes).');
                end
                for i=1:length(elecOri)
                    if strcmpi(elecType{i},'pad')
                        if isempty(elecOri{i})
                            elecOri{i} = 'lr';
                        else
                            if ischar(elecOri{i})
                                if ~any(strcmpi(elecOri{i},{'lr','ap','si'}))
                                    error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                                end
                            else
                                if size(elecOri{i},2)~=3
                                    error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
                                end
                                if size(elecOri{i},1)>1
                                    error('You''re placing more than 1 type of electrodes. Please put orientation info for each pad electrode as a 1-by-3 vector or one of the three keywords ''lr'', ''ap'', or ''si'' in a cell array for option ''elecOri''.');
                                end
                            end
                        end
                    else
                        %                     warning('You''re not placing pad electrodes; customized orientation options will be ignored.');
                        elecOri{i} = [];
                    end
                end
            end
        end
    end
    
    elecPara = struct('capType',capType,'elecType',elecType,...
        'elecSize',elecSize,'elecOri',elecOri);
    
else
    
    fid = fopen('./elec72.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
    elecName = C{4}; for i=1:length(elecName), elecName{i} = strrep(elecName{i},'.',''); end
    capType = '1010';
    elecType = 'disc';
    elecSize = [6 2];
    elecOri = [];
    
    elecPara = struct('capType',capType,'elecType',elecType,...
        'elecSize',elecSize,'elecOri',elecOri);
    
end

if ~exist('T2','var')
    T2 = [];
else
    if ~exist(T2,'file'), error(['The T2 MRI you provided ' T2 ' does not exist.']); end
    
    t2Data = load_untouch_nii(T2);
    if t2Data.hdr.hist.qoffset_x == 0 && t2Data.hdr.hist.srow_x(4)==0
        error('The MRI has a bad header. SPM cannot generate the segmentation properly for MRI with bad header. You can manually align the MRI in SPM Display function to fix the header.');
    end
    % check if bad MRI header    
end

if ~exist('multiaxial','var')
    multiaxial = 0;
else
    if ~ischar(multiaxial), error('Unrecognized option value. Please enter ''on'' or ''off'' for option ''multiaxial''.'); end
    if strcmpi(multiaxial,'off')
        multiaxial = 0;
    elseif strcmpi(multiaxial,'on')
        multiaxial = 1;
    else
        error('Unrecognized option value. Please enter ''on'' or ''off'' for option ''multiaxial''.');
    end
end

if multiaxial && ~isempty(T2)
    error('Multiaxial cannot be run with both T1 and T2 images. If you meant to use Multiaxial, please only provide T1 image with option ''multiaxial'' turned on.');
end

if multiaxial && ~isempty(T2)
    error('Multiaxial cannot be run with both T1 and T2 images. If you meant to use Multiaxial, please only provide T1 image with option ''multiaxial'' turned on.');
end

if ~exist('manual_gui','var')
    manual_gui = 0;
else
    if ~ischar(manual_gui), error('Unrecognized option value. Please enter ''on'' or ''off'' for option ''manual_gui''.'); end
    if strcmpi(manual_gui,'off')
        manual_gui = 0;
    elseif strcmpi(manual_gui,'on')
        manual_gui = 1;
    else
        error('Unrecognized option value. Please enter ''on'' or ''off'' for option ''manual_gui''.');
    end
end

if ~exist('meshOpt','var')
    % meshOpt = struct('radbound',5,'angbound',30,'distbound',0.4,'reratio',3,'maxvol',10);
    meshOpt = struct('radbound',5,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',10);
    % mesh option defaults changed for higher-resolution mesh in version 3
    % meshOpt = struct('radbound',3,'angbound',30,'distbound',0.3,'reratio',3,'maxvol',5);
else
    if ~isstruct(meshOpt), error('Unrecognized format of mesh options. Please enter as a structure, with field names as ''radbound'', ''angbound'', ''distbound'', ''reratio'', and ''maxvol''. Please refer to the iso2mesh documentation for more details.'); end
    meshOptNam = fieldnames(meshOpt);
    if isempty(meshOptNam) || ~all(ismember(meshOptNam,{'radbound';'angbound';'distbound';'reratio';'maxvol'}))
        error('Unrecognized mesh options detected. Supported mesh options are ''radbound'', ''angbound'', ''distbound'', ''reratio'', and ''maxvol''. Please refer to the iso2mesh documentation for more details.');
    end
    if ~isfield(meshOpt,'radbound')
        meshOpt.radbound = 5;
    else
        if ~isnumeric(meshOpt.radbound) || meshOpt.radbound<=0
            error('Please enter a positive number for the mesh option ''radbound''.');
        end
    end
    if ~isfield(meshOpt,'angbound')
        meshOpt.angbound = 30;
    else
        if ~isnumeric(meshOpt.angbound) || meshOpt.angbound<=0
            error('Please enter a positive number for the mesh option ''angbound''.');
        end
    end
    if ~isfield(meshOpt,'distbound')
        meshOpt.distbound = 0.3;
    else
        if ~isnumeric(meshOpt.distbound) || meshOpt.distbound<=0
            error('Please enter a positive number for the mesh option ''distbound''.');
        end
    end
    if ~isfield(meshOpt,'reratio')
        meshOpt.reratio = 3;
    else
        if ~isnumeric(meshOpt.reratio) || meshOpt.reratio<=0
            error('Please enter a positive number for the mesh option ''reratio''.');
        end
    end
    if ~isfield(meshOpt,'maxvol')
        meshOpt.maxvol = 10;
    else
        if ~isnumeric(meshOpt.maxvol) || meshOpt.maxvol<=0
            error('Please enter a positive number for the mesh option ''maxvol''.');
        end
    end
    warning('You''re changing the advanced options of ROAST. Unless you know what you''re doing, please keep mesh options default.');
end

if ~exist('conductivities','var')
    conductivities = struct('white',0.126,'gray',0.276,'csf',1.65,'bone',0.01,...
                           'skin',0.465,'air',2.5e-14,'gel',0.3,'electrode',5.9e7); % literature values
else
    if ~isstruct(conductivities), error('Unrecognized format of conductivity values. Please enter as a structure, with field names as ''white'', ''gray'', ''csf'', ''bone'', ''skin'', ''air'', ''gel'' and ''electrode''.'); end
    conductivitiesNam = fieldnames(conductivities);
    if isempty(conductivitiesNam) || ~all(ismember(conductivitiesNam,{'white';'gray';'csf';'bone';'skin';'air';'gel';'electrode'}))
        error('Unrecognized tissue names detected. Supported tissue names in the conductivity option are ''white'', ''gray'', ''csf'', ''bone'', ''skin'', ''air'', ''gel'' and ''electrode''.');
    end
    if ~isfield(conductivities,'white')
        conductivities.white = 0.126;
    else
        if ~isnumeric(conductivities.white) || any(conductivities.white(:)<=0)
            error('Please enter a positive number for the white matter conductivity.');
        end
        if length(conductivities.white(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
    end
    if ~isfield(conductivities,'gray')
        conductivities.gray = 0.276;
    else
        if ~isnumeric(conductivities.gray) || any(conductivities.gray(:)<=0)
            error('Please enter a positive number for the gray matter conductivity.');
        end
        if length(conductivities.gray(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
    end
    if ~isfield(conductivities,'csf')
        conductivities.csf = 1.65;
    else
        if ~isnumeric(conductivities.csf) || any(conductivities.csf(:)<=0)
            error('Please enter a positive number for the CSF conductivity.');
        end
        if length(conductivities.csf(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
    end
    if ~isfield(conductivities,'bone')
        conductivities.bone = 0.01;
    else
        if ~isnumeric(conductivities.bone) || any(conductivities.bone(:)<=0)
            error('Please enter a positive number for the bone conductivity.');
        end
        if length(conductivities.bone(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
    end
    if ~isfield(conductivities,'skin')
        conductivities.skin = 0.465;
    else
        if ~isnumeric(conductivities.skin) || any(conductivities.skin(:)<=0)
            error('Please enter a positive number for the skin conductivity.');
        end
        if length(conductivities.skin(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
    end
    if ~isfield(conductivities,'air')
        conductivities.air = 2.5e-14;
    else
        if ~isnumeric(conductivities.air) || any(conductivities.air(:)<=0)
            error('Please enter a positive number for the air conductivity.');
        end
        if length(conductivities.air(:))>1, error('Tensor conductivity not supported by ROAST. Please enter a scalar value for conductivity.'); end
    end
    if ~isfield(conductivities,'gel')
        conductivities.gel = 0.3;
    else
        if ~isnumeric(conductivities.gel) || any(conductivities.gel(:)<=0)
            error('Please enter a positive number for the gel conductivity.');
        end
        if length(conductivities.gel(:))>1 && length(conductivities.gel(:))~=length(elecName)
           error('You want to assign different conductivities to the conducting media under different electrodes, but didn''t tell ROAST clearly which conductivity each electrode should use. Please follow the order of electrodes you put in ''recipe'' to give each of them the corresponding conductivity in a vector as the value for the ''gel'' field in option ''conductivities''.');
        end
    end
    if ~isfield(conductivities,'electrode')
        conductivities.electrode = 5.9e7;
    else
        if ~isnumeric(conductivities.electrode) || any(conductivities.electrode(:)<=0)
            error('Please enter a positive number for the electrode conductivity.');
        end
        if length(conductivities.electrode(:))>1 && length(conductivities.electrode(:))~=length(elecName)
           error('You want to assign different conductivities to different electrodes, but didn''t tell ROAST clearly which conductivity each electrode should use. Please follow the order of electrodes you put in ''recipe'' to give each of them the corresponding conductivity in a vector as the value for the ''electrode'' field in option ''conductivities''.');
        end
    end
    warning('You''re changing the advanced options of ROAST. Unless you know what you''re doing, please keep conductivity values default.');
end

if length(conductivities.gel(:))==1
    conductivities.gel = repmat(conductivities.gel,1,length(elecName));
end
if length(conductivities.electrode(:))==1
    conductivities.electrode = repmat(conductivities.electrode,1,length(elecName));
end

if ~exist('simTag','var'), simTag = []; end

if ~exist('doResamp','var')
    doResamp = 0;
else
    if ~ischar(doResamp), error('Unrecognized option value. Please enter ''on'' or ''off'' for option ''resampling''.'); end
    if strcmpi(doResamp,'off')
        doResamp = 0;
    elseif strcmpi(doResamp,'on')
        doResamp = 1;
    else
        error('Unrecognized option value. Please enter ''on'' or ''off'' for option ''resampling''.');
    end
end

if ~exist('paddingAmt','var')
    paddingAmt = 0;
else
    if paddingAmt<=0 || mod(paddingAmt,1)~=0
        error('Unrecognized option value. Please enter positive integer value for option ''zeroPadding''. A recommended value is 10.');
    end
end

% preprocess MRI data
[dirname,subjName,ext] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

if ~strcmpi(subj,'example/nyhead.nii') % only when it's not NY head
    
    t1Data = load_untouch_nii(subj);
    if t1Data.hdr.hist.qoffset_x == 0 && t1Data.hdr.hist.srow_x(4)==0
        error('The MRI has a bad header. SPM cannot generate the segmentation properly for MRI with bad header. You can manually align the MRI in SPM Display function to fix the header.');
    end
    % check if bad MRI header

    if any(t1Data.hdr.dime.pixdim(2:4)<0.8) && ~doResamp
        warning('The MRI has higher resolution (<0.8mm) in at least one direction. This will make the modeling process more computationally expensive and thus slower. If you wish to run faster using just 1-mm model, you can ask ROAST to re-sample the MRI into 1 mm first, by turning on the ''resampling'' option.');
    end
    % check if high-resolution MRI (< 0.8 mm in any direction)
    
    if length(unique(t1Data.hdr.dime.pixdim(2:4)))>1 && ~doResamp
        warning('The MRI has anisotropic resolution. It is highly recommended that you turn on the ''resampling'' option, as the electrode size will not be exact if the model is built from an MRI with anisotropic resolution.');
    end
    % check if anisotropic resolution MRI
    
    [subjRas,isNonRAS] = convertToRAS(subj);
    % check if in non-RAS orientation, and if yes, put it into RAS
    
    [subjRasRS,doResamp] = resampToOneMM(subjRas,doResamp);    
    
    if paddingAmt>0
        subjRasRSPD = zeroPadding(subjRasRS,paddingAmt);
    else
        subjRasRSPD = subjRasRS;
    end
    [~,subjModelName] = fileparts(subjRasRSPD);

    % check for T2, and if yes, check if T2 is aligned with T1
    if ~isempty(T2)
        T2 = realignT2(T2,subjRasRSPD);
        subjRasRSPDspm = [dirname filesep subjModelName '_T1andT2' ext];
    else
        subjRasRSPDspm  = [dirname filesep subjModelName '_T1orT2' ext];
    end
    [~,subjModelNameAftSpm] = fileparts(subjRasRSPDspm);

    % check if multiaxial is on, and if yes, uses it for segmentation
    if multiaxial 
        subjRasRSPDSeg = [dirname filesep subjModelNameAftSpm '_multiaxial' ext];
    else
        %if not, uses SPM for segmentation
        subjRasRSPDSeg = [dirname filesep subjModelNameAftSpm '_SPM' ext];
    end
    [~,subjModelNameAftSeg] = fileparts(subjRasRSPDSeg);

else
    
    if ~exist('example/nyhead_T1orT2_SPM_masks.nii','file')
        unzip('example/nyhead_T1orT2_SPM_masks.nii.zip','example')
    end
    
    isNonRAS = 0; % New York head is in RAS
    
    if doResamp
        error('The beauty of New York head is its 0.5 mm resolution. It''s a bad practice to resample it into 1 mm. Use another head ''example/MNI152_T1_1mm.nii'' for 1 mm model.');
    end
    
    if paddingAmt>0
        zeroPadding('example/nyhead_T1orT2_SPM_masks.nii',paddingAmt);
        subjRasRSPD = ['example/nyhead_padded' num2str(paddingAmt) '.nii'];
        if ~exist(['example/nyhead_padded' num2str(paddingAmt) '_T1orT2_seg8.mat'],'file')
            load('example/nyhead_T1orT2_seg8.mat','image','tpm','Affine');
            origin = inv(image.mat)*[0;0;0;1];
            origin = origin(1:3) + paddingAmt;
            image.mat(1:3,4) = [-dot(origin,image.mat(1,1:3));-dot(origin,image.mat(2,1:3));-dot(origin,image.mat(3,1:3))];
            save(['example/nyhead_padded' num2str(paddingAmt) '_T1orT2_seg8.mat'],'image','tpm','Affine');
        end
    else
        subjRasRSPD = subj;
    end
    [~,subjModelName] = fileparts(subjRasRSPD);

    if ~isempty(T2)
       warning('New York head selected. Any specified T2 image will be ignored.');
       T2 = [];
    end

    if multiaxial
       warning('New York head selected. Multiaxial option will be ignored.');
       multiaxial = 0;
    end

    subjRasRSPDspm = 'example/nyhead_T1orT2.nii';
    subjRasRSPDSeg = 'example/nyhead_T1orT2_SPM.nii';

end

% preprocess electrodes
[elecPara,indInUsrInput] = elecPreproc(subj,elecName,elecPara);

if any(~strcmpi(recipe,'leadfield'))
    
    elecName = elecName(indInUsrInput);
    injectCurrent = injectCurrent(indInUsrInput);
    
    configTxt = [];
    for i=1:length(elecName)
        configTxt = [configTxt elecName{i} ' (' num2str(injectCurrent(i)) ' mA), '];
    end
    configTxt = configTxt(1:end-2);
    
else
    
    elecNameOri = elecName; % back up for re-ordering solutions back to .loc file order;
                            % this is ugly, as .loc file has a different order of electrodes
                            % for historical reasons;
                            % HDE follows .loc file; ROAST follows capInfo.xls
    elecName = elecName(indInUsrInput);
    configTxt = 'leadFieldGeneration';
    
end

conductivities.gel = conductivities.gel(indInUsrInput);
conductivities.electrode = conductivities.electrode(indInUsrInput);

% sort elec options
if length(elecPara)==1
    if size(elecSize,1)>1, elecPara.elecSize = elecPara.elecSize(indInUsrInput,:); end
    if ~ischar(elecOri) && size(elecOri,1)>1
        elecPara.elecOri = elecPara.elecOri(indInUsrInput,:);
    end
elseif length(elecPara)==length(elecName)
    elecPara = elecPara(indInUsrInput);
else
    error('Something is wrong!');
end

options = struct('configTxt',configTxt,'elecPara',elecPara,'T2',T2,'multiaxial',multiaxial,'manual_gui',manual_gui,'meshOpt',meshOpt,'conductivities',conductivities,'uniqueTag',simTag,'resamp',doResamp,'zeroPad',paddingAmt,'isNonRAS',isNonRAS);

% log tracking
Sopt = dir([dirname filesep subjName '_*_roastOptions.mat']);
if isempty(Sopt)
    options = writeRoastLog(subj,options,'roast');
else
    isNew = zeros(length(Sopt),1);
    for i=1:length(Sopt)
        load([dirname filesep Sopt(i).name],'opt');
        isNew(i) = isNewOptions(options,opt,'roast');
    end
    if all(isNew)
        options = writeRoastLog(subj,options,'roast');
    else
        load([dirname filesep Sopt(find(~isNew)).name],'opt');
        if ~isempty(options.uniqueTag) && ~strcmp(options.uniqueTag,opt.uniqueTag)
            warning(['The simulation with the same options has been run before under tag ''' opt.uniqueTag '''. The new tag you specified ''' options.uniqueTag ''' will be ignored.']);
        end
        options.uniqueTag = opt.uniqueTag;
    end
end
uniqueTag = options.uniqueTag;

fprintf('\n');
disp('======================================================')
if ~strcmp(subjName,'nyhead')
    disp(['ROAST ' subj])
else
    disp('ROAST New York head')
end
disp('USING RECIPE:')
disp(configTxt)
disp('...and simulation options saved in:')
disp([dirname filesep subjName '_roastLog,'])
disp(['under tag: ' uniqueTag])
disp('======================================================')
fprintf('\n\n');

% warn users lead field will take a long time to generate
if all(strcmpi(recipe,'leadfield'))
    [~,indRef] = ismember('Iz',elecName);
    indStimElec = setdiff(1:length(elecName),indRef);
    [isInRoastCore,indInRoastCore] = ismember(elecNameOri,elecName(indStimElec));
    isSolved = zeros(length(indStimElec),1);
    for i=1:length(indStimElec)
        if exist([dirname filesep subjName '_' uniqueTag '_e' num2str(indStimElec(i)) '.pos'],'file')
            isSolved(i) = 1;
        end
    end
    % only warn users the first time they run for this subject
    if all(~isSolved) && ~exist([dirname filesep subjName '_' uniqueTag '_roastResult.mat'],'file')
        warning('You specified the ''recipe'' as the ''lead field generation''. Nice choice! Note all customized options on electrodes are overwritten by the defaults. Refer to the readme file for more details. Also this will usually take a long time (>1 day) to generate the lead field for all the candidate electrodes.');
        % doLFconfirm = input('Do you want to continue? ([Y]/N)','s');
        % if strcmpi(doLFconfirm,'n'), disp('Aborted.'); return; end
    end
end

if ~strcmp(subjName,'nyhead')
    if ~multiaxial
        if  ~exist([dirname filesep subjModelNameAftSpm '_seg8.mat'], 'file')
        disp('======================================================')
        disp('        STEP 1 (out of 6): SEGMENT THE MRI...         ')
        disp('======================================================')
        start_seg(subjRasRSPD,T2);
        renameSPMres(subjRasRSPD,subjRasRSPDspm); % rename SPM outputs properly
        else
            disp('======================================================')
            disp('         MRI ALREADY SEGMENTED, SKIP STEP 1           ')
            disp('======================================================')
        end
    end

    if ~exist([dirname filesep subjModelNameAftSeg '_masks.nii'], 'file')
        if multiaxial
            disp('======================================================')
            disp('    STEP 1 (out of 6): MULTIAXIAL SEGMENTATION ...    ')
            disp('======================================================')
            runMultiaxial(subjRasRSPD,subjRasRSPDspm);
        else
            disp('======================================================')
            disp('    STEP 2 (out of 6): SPM SEGMENTATION TOUCHUP ...   ')
            disp('======================================================')
            segTouchup(subjRasRSPDspm,subjRasRSPDSeg);
        end
    else
        if ~multiaxial
        disp('======================================================')
        disp('     SPM SEGMENTATION ALREADY DONE, SKIP STEP 2       ')
        disp('======================================================')
        else 
        disp('======================================================')
        disp('   MULTIAXIAL SEGMENTATION ALREADY DONE, SKIP STEP 1  ')
        disp('======================================================')
        end
    end
else
    disp('======================================================')
    disp(' NEW YORK HEAD SELECTED, GOING TO STEP 3 DIRECTLY...  ')
    disp('======================================================')
    warning('New York head is a 0.5 mm model so is more computationally expensive. Make sure you have a decent machine (>50GB memory) to run ROAST with New York head.')
end
    if ~exist([dirname filesep subjName '_affine_matrix.txt'],"file")
    disp('======================================================')
    disp('    STEP 2 (out of 6): RUNNING NIFTYREG ALIGNMENT ... ')
    disp('======================================================')
    runNiftyReg(subjRas);
    else
    disp('======================================================')
    disp('    NIFTYREG ALIGNMENT ALREADY DONE, SKIP STEP 2 ...  ')
    disp('======================================================')
    end
if ~exist([dirname filesep subjName '_' uniqueTag '_mask_elec.nii'],'file')
    disp('======================================================')
    disp('      STEP 3 (out of 6): ELECTRODE PLACEMENT...       ')
    disp('======================================================')
    hdrInfo = electrodePlacement(subj,subjRasRSPD,subjRasRSPDspm,subjRasRSPDSeg,elecName,options,uniqueTag);
else
    disp('======================================================')
    disp('         ELECTRODE ALREADY PLACED, SKIP STEP 3        ')
    disp('======================================================')
%   load([dirname filesep subjName '_' uniqueTag '_labelVol.mat'],'volume_elecLabel','volume_gelLabel');
    load([dirname filesep subjModelName '_header.mat'],'hdrInfo');
end

if ~exist([dirname filesep subjName '_' uniqueTag '.mat'],'file')
    disp('======================================================')
    disp('        STEP 4 (out of 6): MESH GENERATION...         ')
    disp('======================================================')
    [node,elem,face] = meshByIso2mesh(subj,subjRasRSPDspm,subjRasRSPDSeg,meshOpt,hdrInfo,uniqueTag);
else
    disp('======================================================')
    disp('          MESH ALREADY GENERATED, SKIP STEP 4         ')
    disp('======================================================')
    load([dirname filesep subjName '_' uniqueTag '.mat'],'node','elem','face');
end

if any(~strcmpi(recipe,'leadfield'))

    if ~exist([dirname filesep subjName '_' uniqueTag '_v.pos'],'file')
        disp('======================================================')
        disp('       STEP 5 (out of 6): SOLVING THE MODEL...        ')
        disp('======================================================')
        prepareForGetDP(subj,node,elem,elecName,uniqueTag);
        indElecSolve = 1:length(elecName);
        solveByGetDP(subj,injectCurrent,conductivities,indElecSolve,uniqueTag,'');
    else
        disp('======================================================')
        disp('           MODEL ALREADY SOLVED, SKIP STEP 5          ')
        disp('======================================================')
        %     load([dirname filesep subjName '_' uniqueTag '_elecMeshLabels.mat'],'label_elec');
    end

    if ~exist([dirname filesep subjName '_' uniqueTag '_roastResult.mat'],'file')
        disp('======================================================')
        disp('STEP 6 (final step): SAVING AND VISUALIZING RESULTS...')
        disp('======================================================')
        [vol_all,ef_mag,ef_all] = postGetDP(subj,subjRasRSPDSeg,node,hdrInfo,uniqueTag);
        visualizeRes(subj,subjRasRSPD,subjRasRSPDspm,subjRasRSPDSeg,T2,node,elem,face,injectCurrent,hdrInfo,uniqueTag,0,vol_all,ef_mag,ef_all);
    else
        disp('======================================================')
        disp('  ALL STEPS DONE, LOADING RESULTS FOR VISUALIZATION   ')
        disp('======================================================')
        load([dirname filesep subjName '_' uniqueTag '_roastResult.mat'],'vol_all','ef_mag','ef_all');
        visualizeRes(subj,subjRasRSPD,subjRasRSPDspm,subjRasRSPDSeg,T2,node,elem,face,injectCurrent,hdrInfo,uniqueTag,1,vol_all,ef_mag,ef_all);
    end

else

    if any(~isSolved) && ~exist([dirname filesep subjName '_' uniqueTag '_roastResult.mat'],'file')
        disp('======================================================')
        disp('    STEP 5 (out of 6): GENERATING THE LEAD FIELD...   ')
        disp('           NOTE THIS WILL TAKE SOME TIME...           ')
        disp('======================================================')
        prepareForGetDP(subj,node,elem,elecName,uniqueTag);
        injectCurrent = ones(length(elecName),1); % 1 mA at each candidate electrode
        injectCurrent(indRef) = -1;
        for i=1:length(indStimElec)
            if ~isSolved(i)
                fprintf('\n======================================================\n');
                disp(['SOLVING FOR ELECTRODE ' num2str(i) ' OUT OF ' num2str(length(indStimElec)) ' ...']);
                fprintf('======================================================\n\n');
                indElecSolve = [indStimElec(i) indRef];
                solveByGetDP(subj,injectCurrent,conductivities,indElecSolve,uniqueTag,num2str(indStimElec(i)));
            else
                disp(['ELECTRODE ' num2str(i) ' HAS BEEN SOLVED, SKIPPING...']);
            end
        end
    else
        disp('======================================================')
        disp('       LEAD FIELD ALREADY GENERATED, SKIP STEP 5      ')
        disp('======================================================')
        %     load([dirname filesep subjName '_' uniqueTag '_elecMeshLabels.mat'],'label_elec');
    end

    if ~exist([dirname filesep subjName '_' uniqueTag '_roastResult.mat'],'file')
        disp('========================================================')
        disp('STEP 6 (final step): ASSEMBLING AND SAVING LEAD FIELD...')
        disp('========================================================')
        postGetDP(subj,[],node,hdrInfo,uniqueTag,indStimElec,indInRoastCore(isInRoastCore));
    else
        disp('======================================================')
        disp('         ALL STEPS DONE, READY TO DO TARGETING        ')
        disp(['         FOR SUBJECT ' subj])
        disp(['         USING TAG ' uniqueTag])
        disp('======================================================')
    end   
end

disp('==================ALL DONE ROAST=======================');
