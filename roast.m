function roast(subj,recipe)
% roast(subj,recipe)
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
% The results are saved as "subjName-date-time_result.mat". And you can
% look up the corresponding stimulation config you defined in the log file
% of this subject ("subjName_log"), by using the date-time string.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% November 2017

if nargin<1 || isempty(subj)
    subj = 'example/MNI152_T1_1mm.nii';
end

if nargin<2
    recipe = {'Fp1',1,'P4',-1};
end

if mod(length(recipe),2)~=0
    error('Unrecognized format of your recipe. Please enter as electrodeName-injectedCurrent pair.');
end

elecName = (recipe(1:2:end-1))';
fid = fopen('./BioSemi74.loc'); C = textscan(fid,'%d %f %f %s'); fclose(fid);
elec = C{4};
for i=1:length(elec), elec{i} = strrep(elec{i},'.',''); end
if ~all(ismember(elecName,elec))
    error('Unrecognized electrode names. Please specify electrodes in the 10-10 EEG system only.');
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

[dirname,baseFilename] = fileparts(subj);
if isempty(dirname), dirname = pwd; end

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
    start_seg(subj);
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
    [rnge_elec,rnge_gel] = electrode_placement(subj,'1010',[],[],[],elecName,uniqueTag);
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