function review_roast()
% review_roast()
%
% By Robert Guggenberger at http://www.agricolab.de/
% Nov 2017

[filename, folder] = uigetfile('.mat', 'Pick the one ending on result.mat');
if ~regexp(filename, 'result')
    error('Pick the file ending on result.mat')
    return
end

try
    uniTag   = split(filename,'_');
    uniTag   = uniTag{end-1};
    baseFilename  = split(filename,[uniTag]);
    baseFilename  = baseFilename{1};
    
    load([folder filesep baseFilename uniTag '.mat'])
    load([folder filesep baseFilename uniTag '_rnge.mat'])
    load([folder filesep baseFilename uniTag '_result.mat'])
catch
    error('Could not find all files.')
end
visualizeRes('example/MNI152_T1_1mm.nii',node,elem,face,vol_all,ef_mag,'',uniTag,1);