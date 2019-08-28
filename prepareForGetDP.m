function prepareForGetDP(P,node,elem,elecNeeded,uniTag)
% prepareForGetDP(P,node,elem,elecNeeded,uniTag)
%
% Prepare to solve in getDP
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

% node = node + 0.5; already done right after mesh

% indNode_elecElm = elem(find(elem(:,5) == 8),1:4);
% X = zeros(size(indNode_elecElm,1),3);
% for e = 1:size(indNode_elecElm,1), X(e,:) = mean ( node(indNode_elecElm(e,:),1:3) ); end
% % figure; plot3(X(:,1),X(:,2),X(:,3),'r.');
% 
% Xt = round(X);
% label_elec = volume_elecLabel(sub2ind(size(volume_elecLabel),Xt(:,1),Xt(:,2),Xt(:,3)));
% indBad = find(label_elec==0);
% indGood = find(label_elec>0);
% [~,indOnGood] = map2Points(X(indBad,:),X(indGood,:),'closest'); % using nearest neighbor to fix bad labeling
% label_elec(indBad) = label_elec(indGood(indOnGood));
% 
% indNode_gelElm = elem(find(elem(:,5) == 7),1:4);
% X = zeros(size(indNode_gelElm,1),3);
% for e = 1:size(indNode_gelElm,1), X(e,:) = mean ( node(indNode_gelElm(e,:),1:3) ); end
% % figure; plot3(X(:,1),X(:,2),X(:,3),'r.');
% 
% Xt = round(X);
% label_gel = volume_gelLabel(sub2ind(size(volume_gelLabel),Xt(:,1),Xt(:,2),Xt(:,3)));
% indBad = find(label_gel==0);
% indGood = find(label_gel>0);
% [~,indOnGood] = map2Points(X(indBad,:),X(indGood,:),'closest'); % using nearest neighbor to fix bad labeling
% label_gel(indBad) = label_gel(indGood(indOnGood));
% 
% save([dirname filesep baseFilename '_' uniTag '_elecMeshLabels.mat'],'label_elec','label_gel');

numOfTissue = 6; % hard coded across ROAST.
numOfElec = length(elecNeeded);

element_elecNeeded = cell(numOfElec,1);
area_elecNeeded = zeros(numOfElec,1);

% resolution = mean(hdrInfo.pixdim);
% % mean() here to handle anisotropic rmesolution; ugly. Maybe just
% % resample MRI to isotropic in the very beginning?
% Not needed now when mesh coordinates are in pseudo-world space % ANDY 2019-03-13

warning('off','MATLAB:TriRep:PtsNotInTriWarnId');
for i=1:numOfElec
    
%     if isempty(indNode_elecElm(label_elec==i,:))
%         error(['Electrode ' elecNeeded{i} ' was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.']);
%     end
%         
%     [faces_elec,verts_elec] = freeBoundary(TriRep(indNode_elecElm(label_elec==i,:),node(:,1:3)));
%     [faces_gel,verts_gel] = freeBoundary(TriRep(indNode_gelElm(label_gel==i,:),node(:,1:3)));
    
    indNode_gelElm = elem(find(elem(:,5) == numOfTissue+i),1:4);
    indNode_elecElm = elem(find(elem(:,5) == numOfTissue+numOfElec+i),1:4);
    
    if isempty(indNode_gelElm)
        error(['Gel under electrode ' elecNeeded{i} ' was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.']);
    end
    
    if isempty(indNode_elecElm)
        error(['Electrode ' elecNeeded{i} ' was not meshed properly. Reasons may be: 1) electrode size is too small so the mesher cannot capture it; 2) mesh resolution is not high enough. Consider using bigger electrodes or increasing the mesh resolution by specifying the mesh options.']);
    end
    
    [faces_gel,verts_gel] = freeBoundary(TriRep(indNode_gelElm,node(:,1:3)));
    [faces_elec,verts_elec] = freeBoundary(TriRep(indNode_elecElm,node(:,1:3)));
    
    [~,iE,iG] = intersect(verts_elec,verts_gel,'rows');
    tempTag = ismember(faces_elec,iE);
    % faces_overlap = faces_elec(sum(tempTag,2)==3,:);
    faces_elecOuter = faces_elec(~(sum(tempTag,2)==3),:);
    [~,Loc] = ismember(verts_elec,node(:,1:3),'rows');
    element_elecNeeded{i} = Loc(faces_elecOuter);
    % calculate the surface area
    a = (verts_elec(faces_elecOuter(:, 2),:) - verts_elec(faces_elecOuter(:, 1),:)); %*resolution;
    b = (verts_elec(faces_elecOuter(:, 3),:) - verts_elec(faces_elecOuter(:, 1),:)); %*resolution;
    c = cross(a, b, 2);
    area_elecNeeded(i) = sum(0.5*sqrt(sum(c.^2, 2)));
    
end
if ~exist([dirname filesep baseFilename '_' uniTag '_usedElecArea.mat'],'file')
    save([dirname filesep baseFilename '_' uniTag '_usedElecArea.mat'],'area_elecNeeded');
end

if ~exist([dirname filesep baseFilename '_' uniTag '_ready.msh'],'file')
    
    disp('setting up boundary conditions...');
    
    fid_in = fopen([dirname filesep baseFilename '_' uniTag '.msh']);
    fid_out = fopen([dirname filesep baseFilename '_' uniTag '_ready.msh'],'w');
    
    numOfPart = length(unique(elem(:,5)));
    while ~feof(fid_in)
        s = fgetl(fid_in);
        
        if strcmp(s,'$Elements')
            fprintf(fid_out,'%s\n',s);
            s = fgetl(fid_in);
            numOfElem = str2num(s);
            fprintf(fid_out,'%s\n',num2str(numOfElem+size(cell2mat(element_elecNeeded),1)));
        elseif strcmp(s,'$EndElements')
            ii = 0;
            for j=1:numOfElec
                for i=1:size(element_elecNeeded{j},1)
                    
                    fprintf(fid_out,'%s \n',[num2str(numOfElem+i+ii) ' 2 2 ' num2str(numOfPart+j) ' ' num2str(numOfPart+j) ' ' num2str(element_elecNeeded{j}(i,1)) ' ' num2str(element_elecNeeded{j}(i,2)) ' ' num2str(element_elecNeeded{j}(i,3))]);
                    
                end
                ii = ii + i;
            end
            
            fprintf(fid_out,'%s\n',s);
        else
            fprintf(fid_out,'%s\n',s);
        end
    end
    
    fclose(fid_in);
    fclose(fid_out);
    
end