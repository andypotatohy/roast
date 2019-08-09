function postGetDPforLF(P,node,elem,hdrInfo,indSolved,indInCore,uniTag)
% postGetDPforLF(P,node,elem,hdrInfo,indSolved,indInCore,uniTag)
% 
% Post processing after solving the model for lead field generation. Save
% the results in Matlab format for roast_target to work.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% August 2019

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end
% [~,baseFilenameRSPD] = fileparts(P2);

% node = node + 0.5; already done right after mesh

% convert pseudo-world coordinates back to voxel coordinates for targeting,
% as targeting code works in the voxel space
for i=1:3, node(:,i) = node(:,i)/hdrInfo.pixdim(i); end

indBrain = elem((elem(:,5)==1 | elem(:,5)==2),1:4);
indBrain = unique(indBrain(:));

Atemp = nan(size(node,1),3);
A = nan(length(indBrain),3,length(indSolved));

disp('assembling the lead field...');
for i=1:length(indSolved)
    
    disp(['packing electrode ' num2str(i) ' out of ' num2str(length(indSolved)) ' ...']);
    fid = fopen([dirname filesep baseFilename '_' uniTag '_e' num2str(indSolved(i)) '.pos']);
    fgetl(fid);
    C = textscan(fid,'%d %f %f %f');
    fclose(fid);
    
    Atemp(C{1},:) = cell2mat(C(2:4));
    
    A(:,:,i) = Atemp(indBrain,:);
end

indAdata = find(~isnan(sum(sum(A,3),2))); % make sure no NaN is in matrix A
A = A(indAdata,:,:);

A = reshape(A,length(indBrain)*3,length(indSolved));

locs = node(indBrain,1:3);
locs = locs(indAdata,:); % ...also applies to mesh coordinates

% re-ordering to match the electrode order in .loc file
A = A(:,indInCore);

disp('saving the final results...')
save([dirname filesep baseFilename '_' uniTag '_result.mat'],'A','locs','-v7.3');

disp('======================================================');
disp('The lead field matrix is saved as:');
disp([dirname filesep baseFilename '_' uniTag '_result.mat']);
disp('======================================================');
% disp('You can also find all the results in the following two text files: ');
% disp(['Voltage: ' dirname filesep baseFilename '_' uniTag '_v.pos']);
% disp(['E-field: ' dirname filesep baseFilename '_' uniTag '_e.pos']);
disp('======================================================');
disp('Look up the detailed info for this simulation in the log file: ');
disp([dirname filesep baseFilename '_log']);
disp(['under the simulation tag "' uniTag '".']);
disp('======================================================');
disp('======================================================');
disp('Now you can do targeting by calling: ');
disp('roast_target(subj,simTag,targetCoord,varargin)');
disp('Please refer to the README and roast_target documentation for details.');
disp('======================================================');