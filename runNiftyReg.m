function Affine = NiftyReg(spmOut)
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

    [dirname, mriName, ~] = fileparts(spmOut);
    
    % File paths
    input_mri = spmOut;
    reference_mri = 'example\MNI152_T1_1mm.nii';
    tmp_forward_file = fullfile(dirname, [mriName, '_tmp_forward.txt']);
    final_affine_file = fullfile(dirname, [mriName, '_affine_matrix.txt']);
    niftyreg_path = [pwd '\lib\NiftyReg\'];

    % If affine already exists, just load and return it
    if isfile(final_affine_file)

        try
            Affine = dlmread(final_affine_file);
            if ~isequal(size(Affine), [4, 4])
                error('Matrix is not 4x4.');
            end
            return;
        catch ME
            warning('Failed to read existing affine matrix: %s', ME.message);
            Affine = [];
            return;
        end
    end

    % Set NiftyReg environment
    setenv('PATH', [getenv('PATH') ';' niftyreg_path]);

    % Step 1: Compute forward matrix (temporary)
    cmd = sprintf('"%sreg_aladin" -ref "%s" -flo "%s" -aff "%s"', ...
                  niftyreg_path, reference_mri, input_mri, tmp_forward_file);
    [status, msg] = system(cmd);
    if status ~= 0
        warning('Forward registration failed: %s', msg);
        Affine = [];
        return;
    end

    % Step 2: Invert it and save as final affine
    cmd = sprintf('"%sreg_transform" -invAffine "%s" "%s"', ...
                  niftyreg_path, tmp_forward_file, final_affine_file);
    [status, msg] = system(cmd);
    if status ~= 0
        warning('Inverse affine computation failed: %s', msg);
        Affine = [];
        return;
    end

    % Delete temporary forward file
    if isfile(tmp_forward_file)
        delete(tmp_forward_file);
    end

    % Step 3: Load final affine
    try
        Affine = dlmread(final_affine_file);
        if ~isequal(size(Affine), [4, 4])
            error('Matrix must be 4x4.');
        end
    catch ME
        warning('Failed to read final affine matrix: %s', ME.message);
        Affine = [];
    end
end
