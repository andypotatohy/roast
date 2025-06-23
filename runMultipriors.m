function runMultiaxial(T1,spmOut)
% runMultiaxial(T1,spmOut)
% 
% This calls Multiaxial segmentation that is based on a deep CNN.
% 
% See Hirsch et al 2021 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8209627/)
% for details.
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% April 2024

% warp the TPM to the individual MRI
[dir,t1name] = fileparts(T1);
if isempty(dir), dir = pwd; end
[~,spmOutName] = fileparts(spmOut);

outputFileName = [dir filesep t1name '_indiTPM.nii'];

if ~exist([outputFileName '.gz'],'file')
    disp(['warping TPM to ' T1 ' ...']);
    load([dir filesep spmOutName '_seg8.mat']);

    tpm = spm_load_priors8(tpm);

    d1 = size(tpm.dat{1});
    d1 = d1(1:3);
    M1 = tpm.M;

    Kb = max(lkp);

    d = image.dim(1:3);

    [x1, x2, o] = ndgrid(1:d(1), 1:d(2), 1);
    x3 = 1:d(3);

    prm = [3 3 3 0 0 0];
    Coef = cell(1, 3);
    Coef{1} = spm_bsplinc(Twarp(:,:,:,1), prm);
    Coef{2} = spm_bsplinc(Twarp(:,:,:,2), prm);
    Coef{3} = spm_bsplinc(Twarp(:,:,:,3), prm);

    M = M1\Affine*image.mat;

    QQ = zeros([d(1:3), Kb], 'single');

    for z = 1:length(x3)
        [t1, t2, t3] = defs(Coef, z, MT, prm, x1, x2, x3, M);
        qq = zeros([d(1:2) Kb]);
        b = spm_sample_priors8(tpm, t1, t2, t3);

        for k1 = 1:Kb
            qq(:,:,k1) = b{k1};
        end

        QQ(:,:,z,:) = reshape(qq, [d(1:2), 1, Kb]);
    end

    sQQ = sum(QQ, 4);
    for k1 = 1:Kb
        QQ(:,:,:,k1) = QQ(:,:,:,k1) ./ sQQ;
    end

    hdr = niftiinfo(T1);
    hdr.ImageSize = cat(2, hdr.ImageSize, 6);
    hdr.PixelDimensions = cat(2, hdr.PixelDimensions, 0);
    hdr.Datatype = 'single';
    hdr.BitsPerPixel = 32;
    niftiwrite(QQ, outputFileName, hdr, 'compressed', 1);
else
    disp(['File ' outputFileName '.gz already exists. Warping TPM skipped...']);
end

% Install miniconda and Multiaxial enviornment if it doesn't exist
% And prepare to run Multiaxial python script
str = computer('arch');
switch str
    case 'win64'
        if ~exist('multiaxialEnv', 'dir')
            system([pwd '\lib\multiaxial\setupWindows.bat']);
        else
            disp('Enviornment already exist on Windows. Skipping setup...');
        end

        % Specify the path to libiomp5md.dll (duplicate)
        dllPath = [pwd '\lib\multiaxial\multiaxialENV\Library\bin\libiomp5md.dll'];
        if exist(dllPath, 'file')
            delete(dllPath);
            disp('Duplicate libiomp5md.dll has been deleted');
        end

        pythonExecutable = [pwd '\lib\multiaxial\multiaxialEnv\python.exe'];

    case 'glnxa64'
        setupScript = [pwd '/lib/multiaxial/setupLinux.sh'];

        % Grant execute permissions to the script
        system(['chmod +x ' setupScript]);

        if ~exist('multiaxialEnvLinux', 'dir')
            system([pwd '/lib/multiaxial/setupLinux.sh']);
        else
            disp('Enviornment already exist on Linux. Skipping setup...');
        end

        pythonExecutable = [pwd '/lib/multiaxial/multiaxialEnvLinux/bin/python3'];

    case 'maci64'
        setupScript = [pwd '/lib/multiaxial/setupMac.sh'];

        % Grant execute permissions to the script
        system(['chmod +x ' setupScript]);

        if ~exist('multiaxialEnvMac', 'dir')
            system([pwd '/lib/multiaxial/setupMac.sh']);
        else
            disp('Enviornment already exist on Mac. Skipping setup...');
        end

        pythonExecutable = [pwd '/lib/multiaxial/multiaxialEnvMac/bin/python3'];

    otherwise
        error('Unsupported operating system!');
end

pythonScript = [pwd '/lib/multiaxial/SEGMENT.py'];
configFile = [pwd '/lib/multiaxial/Segmentation_config.py'];

T1 = ['"' T1 '"']; % To allow spaces in subject path

% Construct the command to execute the Python script
command = [pythonExecutable ' ' pythonScript ' ' configFile ' ' T1];

% Execute the Python script
status = system(command);

if status == 0
    disp('Python script executed successfully.');
else
    error('Error running Python script.');
end