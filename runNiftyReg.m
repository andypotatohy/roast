function Affine = runNiftyReg(spmOut)
% Computes and saves the inverse affine matrix as *_affine_matrix.txt using NiftyReg.
% Skips computation if the file already exists.
%
% INPUT:
%   spmOut - Full path to subject MRI NIfTI file
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

[dirname, mriName, ~] = fileparts(spmOut);

% Define paths
input_mri = spmOut;
output_mri = [dirname filesep mriName '_niftyReg.nii'];
reference_mri = fullfile('example', 'MNI152_T1_1mm.nii');
tmp_forward_file = fullfile(dirname, [mriName, '_tmp_forward.txt']);

cmd = sprintf('"%sreg_aladin" -ref "%s" -flo "%s" -aff "%s" -res "%s"', ...
    niftyreg_path, reference_mri, input_mri, tmp_forward_file, output_mri);
[status, msg] = system(cmd);
if status ~= 0
    warning('Forward registration failed: %s', msg);
    Affine = [];
else
    Affine = dlmread(tmp_forward_file);
    delete(tmp_forward_file);
    delete(output_mri);
    Affine = inv(Affine);
    save([dirname filesep mriName '_niftyReg.mat'],'Affine');
end