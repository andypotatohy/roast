function [vol_all,ef_mag] = postGetDP(P1,P2,node,hdrInfo,uniTag)
% [vol_all,ef_mag] = postGetDP(P1,P2,node,hdrInfo,uniTag)
% 
% Post processing after solving the model. Save the result in Matlab format
% in the MRI voxel space.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,baseFilename] = fileparts(P1);
if isempty(dirname), dirname = pwd; end
[~,baseFilenameRSPD] = fileparts(P2);

[xi,yi,zi] = ndgrid(1:hdrInfo.dim(1),1:hdrInfo.dim(2),1:hdrInfo.dim(3));

% node = node + 0.5; already done right after mesh

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
save([dirname filesep baseFilename '_' uniTag '_result.mat'],'vol_all','ef_all','ef_mag');

template = load_untouch_nii(P2);
% Load the original MRI to save the results as NIFTI format
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
disp([dirname filesep baseFilename '_' uniTag '_result.mat']);
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
disp([dirname filesep baseFilename '_log']);
disp(['under the simulation tag "' uniTag '".']);
disp('======================================================');