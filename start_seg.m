function start_seg(P,T2,Template,norm)
% start_seg(P,T2,Template,norm)
%
% Gateway script to start running the SPM8 new segment.
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
% yhuang16@citymail.cuny.edu

if nargin <1 || isempty(P) %no files
    P = spm_select(inf,'image','Select images for new segment');
end

if nargin <2 || isempty(T2) %no T2 specified
    T2 = [];
end

if nargin <3 || isempty(Template) %no Template
    % Template = fullfile(spm('Dir'),'toolbox','Seg','TPM.nii'); %SPM8 default template for new segment
    Template = fullfile(fileparts(which(mfilename)),'eTPM.nii');
    % Template = spm_select(1,'image','Select template for new segment');
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

for i=1:size(P,1)
    
    ref = deblank(P(i,:));
%     [pth,nam,ext] = spm_fileparts(ref);
    [dirname,baseFilename,ext] = fileparts(ref);
    if isempty(dirname), dirname = pwd; end
    
    if strcmp(ext,'.hdr'), ref = [dirname filesep baseFilename '.img']; end
    
    t1Data = load_untouch_nii(ref);
    sliceshow(t1Data.img,[],'gray',[],[],'MRI: Click anywhere to navigate.'); drawnow
    if t1Data.hdr.hist.qoffset_x == 0 && t1Data.hdr.hist.srow_x(4)==0
        error('The MRI has a bad header. SPM cannot generate the segmentation properly for MRI with bad header. You can manually align the MRI in SPM8 Display function to fix the header.');
    end
    
    matlabbatch{1}.spm.tools.preproc8.channel.vols = {ref};  % image to be segmented
    matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001; % P(beta)
    matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
    matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0]; % do NOT save bias field or bias-field corrected image  % ANDY 2013-05-03
    if ~isempty(T2) %T2 specified
        ref2 = deblank(T2(i,:));
%         [pth2,nam2,ext2] = spm_fileparts(ref2);
        [dirname2,baseFilename2,ext2] = fileparts(ref2);
        if isempty(dirname2), dirname2 = pwd; end

        if strcmp(ext2,'.hdr'), ref2 = [dirname2 filesep baseFilename2 '.img']; end
%         fprintf('Using %s to segment T1 and T2 images: %s %s\n', Template, ref, ref2);

        t2Data = load_untouch_nii(ref2);
        sliceshow(t2Data.img,[],'gray',[],[],'MRI: T2. Click anywhere to navigate.'); drawnow
        if t2Data.hdr.hist.qoffset_x == 0 && t2Data.hdr.hist.srow_x(4)==0
            error('The MRI has a bad header. SPM cannot generate the segmentation properly for MRI with bad header. You can manually align the MRI in SPM8 Display function to fix the header.');
        end
        if any(size(t1Data.img)~=size(t2Data.img))
            error('T2 image is not registered to T1 image space.');
        end
        matlabbatch{1}.spm.tools.preproc8.channel(2).vols = {ref2}; % the 2nd image aiding segmentation
        matlabbatch{1}.spm.tools.preproc8.channel(2).biasreg = 0.0001;
        matlabbatch{1}.spm.tools.preproc8.channel(2).biasfwhm = 60;
        matlabbatch{1}.spm.tools.preproc8.channel(2).write = [0 0]; % do NOT save bias field or bias-field corrected image  % ANDY 2013-05-03
    else
%         fprintf('Using %s to segment T1 %s\n', Template, ref);
    end
    
    matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[Tem,',1']}; % TPM tissue #1
    matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2; % number of Gaussians to model this tissue
    matlabbatch{1}.spm.tools.preproc8.tissue(1).native = n; % save the segmented images in native space only ([1 0]) or native+normalized space([1 1])
    matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = w; % save the spatially normalised versions of the segmented images with modulation ([0 1]) or without modulation ([1 0]). We do NOT save anything for this option ([0 0]) at this stage.
    matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[Tem,',2']}; % TPM tissue #2
    matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(2).native = n;
    matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = w;
    matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[Tem,',3']}; % TPM tissue #3
    matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(3).native = n;
    matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = w;
    matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[Tem,',4']}; % TPM tissue #4
    matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.tools.preproc8.tissue(4).native = n;
    matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = w;
    matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[Tem,',5']}; % TPM tissue #5
    matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.tools.preproc8.tissue(5).native = n;
    matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = w;
    matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[Tem, ',6']}; % TPM tissue #6
    matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(6).native = n;
    matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = w;
    
    matlabbatch{1}.spm.tools.preproc8.warp.reg = 4; % P(alpha)
    matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni'; % template type for initial affine registration. Important!
    matlabbatch{1}.spm.tools.preproc8.warp.samp = 3; % downsample distance for segmentation
    matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0]; % save the deformation/warping fields: forward deform ([0 1]), inverse deform ([1 0]), or both ([1 1]). We do NOT save anything for this option ([0 0]) at this stage. % ANDY 2013-05-03
    
    spm_jobman('run',matlabbatch);
    
    if isempty(T2)
        for t=1:6
            movefile([dirname filesep 'c' num2str(t) baseFilename '.nii'],...
                [dirname filesep 'c' num2str(t) baseFilename '_T1orT2.nii']);
        end
        movefile([dirname filesep baseFilename '_rmask.mat'],...
            [dirname filesep baseFilename '_T1orT2_rmask.mat']);
        movefile([dirname filesep baseFilename '_seg8.mat'],...
            [dirname filesep baseFilename '_T1orT2_seg8.mat']);
    else
        for t=1:6
            movefile([dirname filesep 'c' num2str(t) baseFilename '.nii'],...
                [dirname filesep 'c' num2str(t) baseFilename '_T1andT2.nii']);
        end
        movefile([dirname filesep baseFilename '_rmask.mat'],...
            [dirname filesep baseFilename '_T1andT2_rmask.mat']);
        movefile([dirname filesep baseFilename '_seg8.mat'],...
            [dirname filesep baseFilename '_T1andT2_seg8.mat']);
    end
end % for each image
