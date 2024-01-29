function start_seg(P,T2,Template,norm)
% start_seg(P,T2,Template,norm)
%
% Gateway script to start running the SPM12 segment function.
%
% P: input image. If you have a bunch of images to segment, just put the
% image file names vertically using strvcat(), e.g.:
% P = strvcat('/home/user/segProject/img1.nii','/home/user/segProject/img2.nii');
% T2: [optional] another input image, usually a T2 image. To segment multiple images,
% put the file names vertically as in P, make sure each T2 corresponds to
% its corresponding T1 in P.
% If you provide multiple images in P and T2, the program will process them
% one by one, all segmented by the same Template. If you want them to be processed
% by different Template, call this function for each image separately, each time
% using a different Template.
% Template: [optional] tissue probability map to be used (TPM), defaults to extended FOV eTPM.nii
% norm: [optional] will results be saved in normalized (true) or native (false) space?
%
% Adapted from Chris Rorden, 2011.11
% Yu (Andy) Huang, 2013.05
% Yu (Andy) Huang, 2018.09 adapted to SPM12
% yhuang16@citymail.cuny.edu

if nargin <1 || isempty(P) %no files
    P = spm_select(inf,'image','Select images for new segment');
end

if nargin <2 || isempty(T2) %no T2 specified
    T2 = [];
end

if nargin <3 || isempty(Template) %no Template
    Template = fullfile(fileparts(which(mfilename)),'eTPM.nii');
end

if nargin <4 || isempty(norm) %norm not specified... do not normalize data
    norm = false;
end

if norm
    n = [1 1];
    w = [1 0];
else
    n = [1 0];
    w = [0 0];
end

[ptht,namt,extt] = spm_fileparts(deblank(Template(1,:)));
Tem = [ptht,filesep,namt,extt];
spm_jobman('initcfg');

%[dirname, baseFilename, ext] = fileparts(P);
%parts = strsplit(baseFilename, '_');
%newParts = parts(1:end-2);
%resultString = strjoin(newParts, '_');
%Q = [dirname filesep resultString ext];

for i=1:size(P,1)
    
    ref = deblank(P(i,:));
    if strcmp(ref(1),'~') && strcmp(computer('arch'),'glnxa64')
        error('SPM does not recognize ''~'' symbol for Linux home directory. Please provide the full path.');
    end
%     [pth,nam,ext] = spm_fileparts(ref);
    [dirname,baseFilename,ext] = fileparts(ref);
    if isempty(dirname), dirname = pwd; end
    
    if strcmp(ext,'.hdr'), ref = [dirname filesep baseFilename '.img']; end
    
    t1Data = load_untouch_nii(ref);
    sliceshow(t1Data.img,[],'gray',[],[],'MRI: Click anywhere to navigate.'); drawnow
        
    matlabbatch{1}.spm.spatial.preproc.channel.vols = {ref};  % image to be segmented
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; % P(beta) % 0.0001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [0 0]; % do NOT save bias field or bias-field corrected image  % ANDY 2013-05-03
    if ~isempty(T2) %T2 specified
        ref2 = deblank(T2(i,:));
        if strcmp(ref2(1),'~') && strcmp(computer('arch'),'glnxa64')
            error('SPM does not recognize ''~'' symbol for Linux home directory. Please provide the full path.');
        end
%         [pth2,nam2,ext2] = spm_fileparts(ref2);
        [dirname2,baseFilename2,ext2] = fileparts(ref2);
        if isempty(dirname2), dirname2 = pwd; end

        if strcmp(ext2,'.hdr'), ref2 = [dirname2 filesep baseFilename2 '.img']; end
%         fprintf('Using %s to segment T1 and T2 images: %s %s\n', Template, ref, ref2);

        t2Data = load_untouch_nii(ref2);
        sliceshow(t2Data.img,[],'gray',[],[],'MRI: T2. Click anywhere to navigate.'); drawnow
        
        matlabbatch{1}.spm.spatial.preproc.channel(2).vols = {ref2}; % the 2nd image aiding segmentation
        matlabbatch{1}.spm.spatial.preproc.channel(2).biasreg = 0.001; % 0.0001;
        matlabbatch{1}.spm.spatial.preproc.channel(2).biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel(2).write = [0 0]; % do NOT save bias field or bias-field corrected image  % ANDY 2013-05-03
    else
%         fprintf('Using %s to segment T1 %s\n', Template, ref);
    end
    
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {[Tem,',1']}; % TPM tissue #1
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1; %2 % number of Gaussians to model this tissue
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = n; % save the segmented images in native space only ([1 0]) or native+normalized space([1 1])
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = w; % save the spatially normalised versions of the segmented images with modulation ([0 1]) or without modulation ([1 0]). We do NOT save anything for this option ([0 0]) at this stage.
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {[Tem,',2']}; % TPM tissue #2
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1; %2
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = n;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = w;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {[Tem,',3']}; % TPM tissue #3
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = n;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = w;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {[Tem,',4']}; % TPM tissue #4
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = n;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = w;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {[Tem,',5']}; % TPM tissue #5
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = n;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = w;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {[Tem, ',6']}; % TPM tissue #6
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = n;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = w;
    
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2]; % P(alpha) % 4;
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni'; % template type for initial affine registration. Important!
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 3; % downsample distance for segmentation
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0]; % save the deformation/warping fields: forward deform ([0 1]), inverse deform ([1 0]), or both ([1 1]). We do NOT save anything for this option ([0 0]) at this stage. % ANDY 2013-05-03

    % below new options from SPM12
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 0; % in spm8, no MRF cleanup was done
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 0; % do not use cleanup, use my own cleanup routine
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    
    spm_jobman('run',matlabbatch);
    
    %[dir, baseFile, ~] = fileparts(P);
    %if isempty(T2)
    %    for t=1:6
    %        movefile([dirname filesep 'c' num2str(t) baseFilename '.nii'],...
    %            [dir filesep 'c' num2str(t) baseFile '.nii']);
    %    end
    %    movefile([dirname filesep baseFilename '_rmask.mat'],...
    %        [dir filesep baseFile '_rmask.mat']);
    %    movefile([dirname filesep baseFilename '_seg8.mat'],...
    %        [dir filesep baseFile '_seg8.mat']);
   % else
    %    for t=1:6
    %        movefile([dirname filesep 'c' num2str(t) baseFilename '.nii'],...
     %           [dir filesep 'c' num2str(t) baseFile '.nii']);
    %    end
     %   movefile([dirname filesep baseFilename '_rmask.mat'],...
     %       [dir filesep baseFile '_rmask.mat']);
     %   movefile([dirname filesep baseFilename '_seg8.mat'],...
     %       [dir filesep baseFile '_seg8.mat']);
    %end
end % for each image