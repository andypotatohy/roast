function runNiftyReg(input_mri)
% Computes and saves the inverse affine matrix as *_affine_matrix.txt using NiftyReg.
% Skips computation if the file already exists.
%
% INPUT:
%   input_mri - Full path to subject MRI NIfTI file
%
% OUTPUT:
%   Affine - 4x4 inverse affine matrix (double)
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
% June 2025

% Choose correct NiftyReg path based on OS
arch = computer('arch');
switch arch
    case 'win64'
        niftyreg_path = fullfile(pwd, 'lib', 'NiftyReg', 'win', filesep);
    case 'glnxa64'
        niftyreg_path = fullfile(pwd, 'lib', 'NiftyReg', 'linux', filesep);
    case 'maci64'
        niftyreg_path = fullfile(pwd, 'lib', 'NiftyReg', 'mac', filesep);
    otherwise
        error('Unsupported operating system!');
end

[dirname, mriName] = fileparts(input_mri);

% Define paths
output_mri = [dirname filesep mriName '_niftyReg.nii'];
reference_mri = fullfile('example', 'MNI152_T1_1mm.nii');
tmp_forward_file = fullfile(dirname, [mriName, '_tmp_forward.txt']);

cmd = sprintf('"%sreg_aladin" -ref "%s" -flo "%s" -aff "%s" -res "%s" -voff', ...
    niftyreg_path, reference_mri, input_mri, tmp_forward_file, output_mri);
disp('This may take a couple of minutes ...');
status = system(cmd);
if status ~= 0
    error('niftyReg failed');
else
    Affine = dlmread(tmp_forward_file);
    delete(tmp_forward_file);
    delete(output_mri);
    Affine = inv(Affine); % to be consistent with SPM format
    % also save voxel-to-world matrices for the MRI and TPM, to be
    % consistent with SPM format
    tpm = spm_vol('eTPM.nii');
    image = spm_vol(input_mri);
    save([dirname filesep mriName '_niftyReg.mat'],'Affine','image','tpm');
end