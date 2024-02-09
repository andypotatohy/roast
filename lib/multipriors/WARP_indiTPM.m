function WARP_indiTPM(subj)

    % Split the path into directory, filename, and extension
    [pth, name, ext] = fileparts(subj);
    pth = strrep(pth, '\', '/');
    % Check for the existence of the output file
    outputFileName = [pth '/' name '_indiTPM.nii.gz'];
    
    % Replace backslashes with forward slashes
    
    if isfile(outputFileName) == 0
        disp(['warping TPM to ' pth '/' name ext ' ...']);
        load([pth '/' name '_seg8.mat']);
    
        % try
        %     image = image(1);
        % catch
        %     image = image;
        % end
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

        hdr = niftiinfo([pth filesep name ext]);
        hdr.ImageSize = cat(2, hdr.ImageSize, 6);
        hdr.PixelDimensions = cat(2, hdr.PixelDimensions, 0);
        hdr.Datatype = 'single';
        hdr.BitsPerPixel = 32;
        niftiwrite(QQ, [pth filesep name '_indiTPM.nii'], hdr, 'compressed', 1);

    else
        disp(['File ' outputFileName ' already exists. Warping TPM skipped...']);
    end
end