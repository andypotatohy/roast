function [para,ind2UI] = elecPreproc(subj,elec,para)
% [para,ind2UI] = elecPreproc(subj,elec,para)
% 
% Pre-process electrode names.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

doPredefined = 0;
doNeck = 0;
doCustom = 0;
unknownElec = 0;

switch lower(para(1).capType)
    case {'1020','1010','1005'}
        load('./cap1005FullWithExtra.mat','capInfo');
        elecPool_P = capInfo{1};
    case 'biosemi'
        load('./capBioSemiFullWithExtra.mat','capInfo');
        elecPool_P = capInfo{1};
end

elecPool_N = {'nk1';'nk2';'nk3';'nk4'};

for i=1:length(elec)
    if ismember(elec{i},elecPool_P)
        doPredefined = 1;
    elseif ismember(lower(elec{i}),elecPool_N)
        doNeck = 1;
    elseif ~isempty(strfind(lower(elec{i}),'custom'))
        doCustom = 1;
        if ~exist('elecPool_C','var')
            [dirname,baseFilename] = fileparts(subj);
            if isempty(dirname), dirname = pwd; end
            fid = fopen([dirname filesep baseFilename '_customLocations']);
            if fid==-1
                error('You specified customized electrode locations but did not provide the location coordinates. Please put together all the coordinates in a text file ''subjectName_customLocations'' and store it under the subject folder.');
            end
            capInfo_C = textscan(fid,'%s %f %f %f');
            elecPool_C = capInfo_C{1};
            fclose(fid);
        end
        if ~ismember(elec{i},elecPool_C)
            fprintf('Unrecognized electrode %s.\n',elec{i});
            unknownElec = unknownElec+1;
        end
    else
        fprintf('Unrecognized electrode %s.\n',elec{i});
        unknownElec = unknownElec+1;
    end
end

if unknownElec>0
    error('Unrecognized electrodes found. It may come from the following mistakes: 1) you specified one cap type (e.g. 1010) but asked the electrode name in the other system (e.g. BioSemi); 2) you defined some customized electrode location but forgot to put ''custom'' as a prefix in the electrode name; 3) you picked up one of the electrodes that falls on the ears or eyes (which are removed, see capInfo.xls); 4) you asked ROAST to do an electrode that does not belong to any system (neither 1005, BioSemi, nor your customized electrodes).');
end

if doPredefined
    [isPredefined,indPredefined] = ismember(elec,elecPool_P);
    ind2UI_P = find(isPredefined);
    indP = indPredefined(isPredefined);
    [indP,indtemp] = sort(indP);
    ind2UI_P = ind2UI_P(indtemp);
else
    indP = [];
    ind2UI_P = [];
end

if doNeck
    [isNeck,indNeck]=ismember(lower(elec),elecPool_N);
    ind2UI_N = find(isNeck);
    indN = indNeck(isNeck);
    [indN,indtemp] = sort(indN);
    ind2UI_N = ind2UI_N(indtemp);
else
    indN = [];
    ind2UI_N = [];
end

if doCustom
    [isCustom,indCustom]=ismember(elec,elecPool_C);
    ind2UI_C = find(isCustom);
    indC = indCustom(isCustom);
    [indC,indtemp] = sort(indC);
    ind2UI_C = ind2UI_C(indtemp);
else
    indC = [];
    ind2UI_C = [];
end

for i=1:length(para)
    para(i).indP = indP; para(i).indN = indN; para(i).indC = indC;
end
if isempty(indP) && isempty(indN) && ~isempty(indC)
    for i=1:length(para), para(i).capType = 'none'; end
end

ind2UI = cat(1,ind2UI_P,ind2UI_N,ind2UI_C);