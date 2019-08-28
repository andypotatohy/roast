function [vol_all,ef_mag,ef_all] = postGetDP(P1,P2,node,hdrInfo,uniTag,indSolved,indInCore)
% [vol_all,ef_mag,ef_all] = postGetDP(P1,P2,node,hdrInfo,uniTag,indSolved,indInCore)
%
% Post processing after solving the model / generating the lead field.
% Save the result in Matlab format in the MRI voxel space. For the lead
% field, it's saved in Matlab format for roast_target() to work.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017
% August 2019 adding lead field

[dirname,baseFilename] = fileparts(P1);
if isempty(dirname), dirname = pwd; end

% node = node + 0.5; already done right after mesh

if ~isempty(P2) % for roast()
    
    % convert pseudo-world coordinates back to voxel coordinates for
    % interpolation into regular grid in the voxel space
    for i=1:3, node(:,i) = node(:,i)/hdrInfo.pixdim(i); end

    [~,baseFilenameRSPD] = fileparts(P2);
    
    [xi,yi,zi] = ndgrid(1:hdrInfo.dim(1),1:hdrInfo.dim(2),1:hdrInfo.dim(3));
    
    disp('converting the results into Matlab format...');
    fid = fopen([dirname filesep baseFilename '_' uniTag '_v.pos']);
    fgetl(fid);
    C = textscan(fid,'%d %f');
    fclose(fid);
    
    C{2} = C{2} - min(C{2}); % re-reference the voltage
    
    F = TriScatteredInterp(node(C{1},1:3), C{2});
    vol_all = F(xi,yi,zi);
    
    fid = fopen([dirname filesep baseFilename '_' uniTag '_e.pos']);
    fgetl(fid);
    C = textscan(fid,'%d %f %f %f');
    fclose(fid);
    
    ef_all = zeros([hdrInfo.dim 3]);
    F = TriScatteredInterp(node(C{1},1:3), C{2});
    ef_all(:,:,:,1) = F(xi,yi,zi);
    F = TriScatteredInterp(node(C{1},1:3), C{3});
    ef_all(:,:,:,2) = F(xi,yi,zi);
    F = TriScatteredInterp(node(C{1},1:3), C{4});
    ef_all(:,:,:,3) = F(xi,yi,zi);
    
    ef_mag = sqrt(sum(ef_all.^2,4));
    
    disp('saving the final results...')
    save([dirname filesep baseFilename '_' uniTag '_roastResult.mat'],'vol_all','ef_all','ef_mag','-v7.3');
    
    if isempty(strfind(P2,'example/nyhead'))
        template = load_untouch_nii(P2);
    else
        template = load_untouch_nii([dirname filesep baseFilenameRSPD '_T1orT2_masks.nii']);
    end % Load the original MRI to save the results as NIFTI format
    
    template.hdr.dime.datatype = 16;
    template.hdr.dime.bitpix = 32;
    template.hdr.dime.scl_slope = 1; % so that display of NIFTI will not alter the data
    template.hdr.dime.cal_max = 0;
    template.hdr.dime.cal_min = 0;
    
    template.img = single(vol_all);
    template.hdr.dime.glmax = max(vol_all(:));
    template.hdr.dime.glmin = min(vol_all(:));
    template.hdr.hist.descrip = 'voltage';
    template.fileprefix = [dirname filesep baseFilename '_' uniTag '_v'];
    save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_v.nii']);
    
    template.img = single(ef_mag);
    template.hdr.dime.glmax = max(ef_mag(:));
    template.hdr.dime.glmin = min(ef_mag(:));
    template.hdr.hist.descrip = 'EF mag';
    template.fileprefix = [dirname filesep baseFilename '_' uniTag '_emag'];
    save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_emag.nii']);
    
    template.hdr.dime.dim(1) = 4;
    template.hdr.dime.dim(5) = 3;
    template.img = single(ef_all);
    template.hdr.dime.glmax = max(ef_all(:));
    template.hdr.dime.glmin = min(ef_all(:));
    template.hdr.hist.descrip = 'EF';
    template.fileprefix = [dirname filesep baseFilename '_' uniTag '_e'];
    save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_e.nii']);
    
    disp('======================================================');
    disp('Results are saved as:');
    disp([dirname filesep baseFilename '_' uniTag '_roastResult.mat']);
    disp('...and also saved as NIFTI files:');
    disp(['Voltage: ' dirname filesep baseFilename '_' uniTag '_v.nii']);
    disp(['E-field: ' dirname filesep baseFilename '_' uniTag '_e.nii']);
    disp(['E-field magnitude: ' dirname filesep baseFilename '_' uniTag '_emag.nii']);
    disp('======================================================');
    disp('You can also find all the results in the following two text files: ');
    disp(['Voltage: ' dirname filesep baseFilename '_' uniTag '_v.pos']);
    disp(['E-field: ' dirname filesep baseFilename '_' uniTag '_e.pos']);
    disp('======================================================');
    disp('Look up the detailed info for this simulation in the log file: ');
    disp([dirname filesep baseFilename '_roastLog']);
    disp(['under the simulation tag "' uniTag '".']);
    disp('======================================================');
    
else % for roast_target()
    
%     indBrain = elem((elem(:,5)==1 | elem(:,5)==2),1:4);
%     indBrain = unique(indBrain(:));
%     
%     Atemp = nan(size(node,1),3);
%     A = nan(length(indBrain),3,length(indSolved));

    A_all = nan(size(node,1),3,length(indSolved));
    
    disp('assembling the lead field...');
    for i=1:length(indSolved)
        
        disp(['packing electrode ' num2str(i) ' out of ' num2str(length(indSolved)) ' ...']);
        fid = fopen([dirname filesep baseFilename '_' uniTag '_e' num2str(indSolved(i)) '.pos']);
        fgetl(fid);
        C = textscan(fid,'%d %f %f %f');
        fclose(fid);
        
%         Atemp(C{1},:) = cell2mat(C(2:4));
%         
%         A(:,:,i) = Atemp(indBrain,:);
        
        A_all(C{1},:,i) = cell2mat(C(2:4));
        
        % to save disk space
        delete([dirname filesep baseFilename '_' uniTag '_e' num2str(indSolved(i)) '.pos']);
    end
    
%     indAdata = find(~isnan(sum(sum(A,3),2))); % make sure no NaN is in matrix A
%     A = A(indAdata,:,:);
%     
%     A = reshape(A,length(indBrain)*3,length(indSolved)); % this is bug
%     
%     locs = node(indBrain,1:3);
%     locs = locs(indAdata,:); % ...also applies to mesh coordinates
    
    % re-ordering to match the electrode order in .loc file
%     A = A(:,indInCore);
    A_all = A_all(:,:,indInCore);
    
    disp('saving the final results...')
%     save([dirname filesep baseFilename '_' uniTag '_roastResult.mat'],'A','locs','-v7.3');
    save([dirname filesep baseFilename '_' uniTag '_roastResult.mat'],'A_all','-v7.3');
    
    disp('======================================================');
    disp('The lead field matrix is saved as:');
    disp([dirname filesep baseFilename '_' uniTag '_roastResult.mat']);
    disp('======================================================');
    % disp('You can also find all the results in the following two text files: ');
    % disp(['Voltage: ' dirname filesep baseFilename '_' uniTag '_v.pos']);
    % disp(['E-field: ' dirname filesep baseFilename '_' uniTag '_e.pos']);
    disp('======================================================');
    disp('Look up the detailed info for this simulation in the log file: ');
    disp([dirname filesep baseFilename '_roastLog']);
    disp(['under the simulation tag "' uniTag '".']);
    disp('======================================================');
    disp('======================================================');
    disp('Now you can do targeting by calling: ');
    disp('roast_target(subj,simTag,targetCoord,varargin)');
    disp('Please refer to the README and roast_target() documentation for details.');
    disp('======================================================');
    
end