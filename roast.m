function roast(subj,recipe,varargin)
% roast(subj,recipe,varargin)
%
% Main function for ROAST.
%
% Example: roast
%
% Default call of ROAST, will demo a modeling process on the MNI152 head.
% Specifically, this will use the MRI of the MNI152 head to build a model
% of transcranial electric stimulation (TES) with anode on Fp1 (1 mA) and cathode
% on P4 (-1 mA).
%
% Example: roast('example/subject1.nii')
%
% Build a TES model on any subject you want, just provide the MRI in the
% 1st argument. The default stimulation config is anode on Fp1 (1 mA) and
% cathode on P4 (-1 mA).
%
% Example: roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4})
%
% Build the TES model on any subject with your own "recipe". Here we inject
% 0.3 mA at electrode F1, 0.7 mA at P2, and we ask 0.6 mA coming out of C5,
% and 0.4 mA flowing out of O2. You can define any stimulation montage you want
% in the 2nd argument, with electrodeName-injectedCurrent pair. The electrodes
% supported in this version come from the 10-10 EEG system, and you can find
% their names in the file 1010electrodes.png under the root directory of ROAST.
% Note the unit of the injected current is milliampere, and make sure they sum
% up to 0.
%
% OPTIONS HERE.
%
% The results are saved as "subjName-date-time_result.mat". And you can
% look up the corresponding stimulation config you defined in the log file
% of this subject ("subjName_log"), by using the date-time string.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% March 2018

if nargin<1 || isempty(subj)
    subj = 'example/MNI152_T1_1mm.nii';
end

if nargin<2
    recipe = {'Fp1',1,'P4',-1};
end

if mod(length(recipe),2)~=0
    error('Unrecognized format of your recipe. Please enter as electrodeName-injectedCurrent pair.');
end

lenOptions = length(varargin);
% if mod(lenOptions,2)~=0
%     error('Unrecognized format of options. Please enter as property-value pair.');
% end

% take in user-specified options
% for i=1:2:lenOptions-1
indArg = 1;
while indArg <=lenOptions
    switch varargin{indArg}
        case 'capType'
            capType = varargin{indArg+1};
            indArg = indArg+2;
        case 'elecType'
            elecType = varargin{indArg+1};
            indArg = indArg+2;
        case 'elecSize'
            elecSize = varargin{indArg+1};
            indArg = indArg+2;
        case 'elecOri'
            elecOri = varargin{indArg+1};
            indArg = indArg+2;
        case 'legacy'
            legacy = 1;
            indArg = indArg+1;
        case 'T2'
            T2 = varargin{indArg+1};
            indArg = indArg+2;
        case 'meshOptions'
            meshOpt = varargin{indArg+1};
            indArg = indArg+2;
        otherwise
            error('Supported options are: ''capType'', ''elecType'', ''elecSize'', ''elecOri'', ''legacy'', ''T2'' and ''meshOptions''.');
    end
end

% set up defaults
if ~exist('capType','var'), capType = '1010'; end
if ~exist('elecType','var'), elecType = 'disc'; end
if ~exist('elecSize','var'), elecSize = [6 2]; end
if ~exist('elecOri','var') && strcmp(elecType,'pad')
    elecOri = 'lr';
else
    elecOri = [];
end
if ~exist('legacy','var'), legacy = 0; end
if ~exist('T2','var'), T2 = []; end
if ~exist('meshOpt','var')
    meshOpt = struct('radbound',5,'angbound',30,'distbound',0.4,'reratio',3,'maxvol',10);
end

% more input judge on option conflicts here
if ischar(elecOri) && ~strcmp(elecOri,'lr') && ~strcmp(elecOri,'ap') && ~strcmp(elecOri,'si')
    error('Unrecognized pad orientation. Please enter ''lr'', ''ap'', or ''si'' for pad orientation; or just enter the direction vector of the long axis of the pad');
end


[dirname,baseFilename] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

switch capType
    case {'1020','1010','1005'}
        load('./cap1005FullWithExtra.mat','capInfo');
        elecPool = capInfo{1};
    case 'biosemi'
        load('./capBioSemiFullWithExtra.mat','capInfo');
        elecPool = capInfo{1};
        %     case 'customized'
        %         fid = fopen([dirname filesep 'customLocations']);
        %         if fid==-1
        %             error('You specified customized electrode locations but did not provide the location coordinates. Please put together all the coordinates in a file ''customLocations'' and store it under the subject folder');
        %         end
        %         capInfo = textscan(fid,'%s %f %f %f');
        %         fclose(fid);
        %         elecPool = capInfo{1};
    otherwise
        %         error('Supported cap types are: ''1020'', ''1010'', ''1005'', ''BioSemi'',
        %         and any locations you specified in file ''customLocations'' YOU stored under the subject folder.');
end
% elecPool = cat(1,elecPool,{'Nk1';'Nk2';'Nk3';'Nk4'});

elecName = (recipe(1:2:end-1))';
% fid = fopen('./BioSemi74.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
% elec = C{4};
% for i=1:length(elec), elec{i} = strrep(elec{i},'.',''); end
% if ~all(ismember(elecName,elec))
%     error('Unrecognized electrode names. Please specify electrodes in the 10-10 EEG system only.');
% end

doPredefined = 0;
doNeck = 0;
doCustom = 0;
unknownElec = 0;
for i=1:length(elecName)
    if ismember(elecName{i},elecPool)
        doPredefined = 1;
    elseif ismember(elecName{i},{'Nk1';'Nk2';'Nk3';'Nk4'})
        doNeck = 1;
    elseif ~isempty(strfind(elecName{i},'custom')) || ~isempty(strfind(elecName{i},'Custom'))
        doCustom = 1;
    else
        fprintf('Unrecognized electrode %s.\n',elecName{i});
        unknownElec = unknownElec+1;
    end
end

% foundElec = ismember(elecName,elecPool);
% indCustomElec = find(~foundElec);
% isUnknownElec = zeros(length(indCustomElec),1);
% for i=1:length(indCustomElec)
%     if isempty(strfind(elecName{indCustomElec(i)},'custom')) && isempty(strfind(elecName{indCustomElec(i)},'Custom'))
%         fprintf('Unrecognized electrode %s.\n',elecName{indCustomElec(i)});
%         isUnknownElec(i) = 1;
%     end
% end
% if any(isUnknownElec)
%     error('Unrecognized electrodes found. It''s possible that you specified one cap type (e.g. 1010) but provided the electrode name in the other system (e.g. BioSemi); it''s also possible that you defined some customized electrode location but forgot to put ''custom'' as a prefix in the electrode name.');
% end
if unknownElec>0
    error('Unrecognized electrodes found. It''s possible that you specified one cap type (e.g. 1010) but provided the electrode name in the other system (e.g. BioSemi); it''s also possible that you defined some customized electrode location but forgot to put ''custom'' as a prefix in the electrode name. If you specified neck electrode, please write ''Nk#'' instead of ''nk#'' or ''NK#''.');
end
if doCustom
    fid = fopen([dirname filesep baseFilename '_customLocations']);
    if fid==-1
        error('You specified customized electrode locations but did not provide the location coordinates. Please put together all the coordinates in a text file ''subjectName_customLocations'' and store it under the subject folder');
    end
    fclose(fid);
end

injectCurrent = cell2mat(recipe(2:2:end));
if sum(injectCurrent)~=0
    error('Electric currents going in and out of the head not balanced. Please make sure they sum to 0.');
end

configTxt = [];
for i=1:length(elecName)
    configTxt = [configTxt elecName{i} ' (' num2str(injectCurrent(i)) ' mA), '];
end
configTxt = configTxt(1:end-2);

% save simulation options in a text file named by uniqueTag


if ~exist([dirname filesep baseFilename '_log'],'file')
    fid = fopen([dirname filesep baseFilename '_log'],'w');
    uniqueTag = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss'));
    % "datetime" only available from Matlab 2014b
    fprintf(fid,'%s\t%s\n',uniqueTag,configTxt);
    fclose(fid);
else
    fid = fopen([dirname filesep baseFilename '_log'],'r');
    C = textscan(fid,'%s\t%s','delimiter','\t');
    fclose(fid);
    [Lia,Loc] = ismember(configTxt,C{2});
    % should use option text file to judge if it's run before
    if ~Lia
        fid = fopen([dirname filesep baseFilename '_log'],'a');
        uniqueTag = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss'));
        % "datetime" only available from Matlab 2014b
        fprintf(fid,'%s\t%s\n',uniqueTag,configTxt);
        fclose(fid);
    else
        uniqueTag = C{1}{Loc};
    end
end

fprintf('\n\n');
disp('======================================================')
disp(['ROAST ' subj])
disp('USING RECIPE:')
disp(configTxt)
disp('======================================================')
fprintf('\n\n');

addpath(genpath([fileparts(which(mfilename)) filesep 'lib/']));

if ~exist([dirname filesep 'c1' baseFilename '.nii'],'file')
    disp('======================================================')
    disp('            STEP 1: SEGMENT THE MRI...                ')
    disp('======================================================')
    start_seg(subj,T2);
else
    disp('======================================================')
    disp('          MRI ALREADY SEGMENTED, SKIP STEP 1          ')
    disp('======================================================')
end

if ~exist([dirname filesep baseFilename '_mask_gray.nii'],'file')
    disp('======================================================')
    disp('          STEP 2: SEGMENTATION TOUCHUP...             ')
    disp('======================================================')
    mysegment(subj);
    autoPatching(subj);
else
    disp('======================================================')
    disp('    SEGMENTATION TOUCHUP ALREADY DONE, SKIP STEP 2    ')
    disp('======================================================')
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '_rnge.mat'],'file')
    disp('======================================================')
    disp('          STEP 3: ELECTRODE PLACEMENT...              ')
    disp('======================================================')
    elecPara = struct('capType',capType,'elecType',elecType,...
                      'elecSize',elecSize,'elecOri',elecOri,'legacy',legacy,...
                      'doPredefined',doPredefined,'doNeck',doNeck,'doCustom',doCustom);
    [rnge_elec,rnge_gel] = electrodePlacement(subj,elecName,elecPara,uniqueTag);
else
    disp('======================================================')
    disp('         ELECTRODE ALREADY PLACED, SKIP STEP 3        ')
    disp('======================================================')
    load([dirname filesep baseFilename '_' uniqueTag '_rnge.mat'],'rnge_elec','rnge_gel');
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '.mat'],'file')
    disp('======================================================')
    disp('            STEP 4: MESH GENERATION...                ')
    disp('======================================================')
    [node,elem,face] = meshByIso2mesh(subj,uniqueTag);
else
    disp('======================================================')
    disp('          MESH ALREADY GENERATED, SKIP STEP 4         ')
    disp('======================================================')
    load([dirname filesep baseFilename '_' uniqueTag '.mat'],'node','elem','face');
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '_v.pos'],'file')
    disp('======================================================')
    disp('           STEP 5: SOLVING THE MODEL...               ')
    disp('======================================================')
    prepareForGetDP(subj,node,elem,rnge_elec,rnge_gel,elecName,uniqueTag);
    solveByGetDP(subj,injectCurrent,uniqueTag);
else
    disp('======================================================')
    disp('           MODEL ALREADY SOLVED, SKIP STEP 5          ')
    disp('======================================================')
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '_result.mat'],'file')
    disp('======================================================')
    disp('     STEP 6: SAVING AND VISUALIZING THE RESULTS...    ')
    disp('======================================================')
    [vol_all,ef_mag] = postGetDP(subj,node,uniqueTag);
    visualizeRes(subj,node,elem,face,vol_all,ef_mag,injectCurrent,uniqueTag,0);
else
    disp('======================================================')
    disp('  ALL STEPS DONE, LOADING RESULTS FOR VISUALIZATION   ')
    disp('======================================================')
    load([dirname filesep baseFilename '_' uniqueTag '_result.mat'],'vol_all','ef_mag');
    visualizeRes(subj,node,elem,face,vol_all,ef_mag,injectCurrent,uniqueTag,1);
end