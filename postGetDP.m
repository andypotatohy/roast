function [vol_all,ef_mag] = postGetDP(P,node,uniTag)
% [vol_all,ef_mag] = postGetDP(P,node,uniTag)
% 
% Post processing after solving the model. Save the result in Matlab format
% in the MRI voxel space.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

if ~strcmp(baseFilename,'nyhead')
    data = load_untouch_nii(P);
else
    data = load_untouch_nii([dirname filesep baseFilename '_T1orT2_mask_skin.nii']);
end
[xi,yi,zi] = ndgrid(1:size(data.img,1), 1:size(data.img,2), 1:size(data.img,3));

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

ef_all = zeros([size(data.img) 3]);
F = TriScatteredInterp(node(C{1},1:3), C{2});
ef_all(:,:,:,1) = F(xi,yi,zi);
F = TriScatteredInterp(node(C{1},1:3), C{3});
ef_all(:,:,:,2) = F(xi,yi,zi);
F = TriScatteredInterp(node(C{1},1:3), C{4});
ef_all(:,:,:,3) = F(xi,yi,zi);

ef_mag = sqrt(sum(ef_all.^2,4));

disp('saving the final results...')
save([dirname filesep baseFilename '_' uniTag '_result.mat'],'vol_all','ef_all','ef_mag');
disp('======================================================');
disp(['Results saved as ' dirname filesep baseFilename '_' uniTag '_result.mat'])
disp('You can also find the results in the two text files: ');
disp([dirname filesep baseFilename '_' uniTag '_v.pos']);
disp(['and ' dirname filesep baseFilename '_' uniTag '_e.pos']);
disp('Look up the detailed info for this simulation in the log file');
disp('by using the unique simulation tag.');
disp('======================================================');