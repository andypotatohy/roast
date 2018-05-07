function [label_elec,label_gel] = prepareForGetDP(P,node,elem,rnge_elec,rnge_gel,hdrInfo,elecNeeded,uniTag)
% [label_elec,label_gel] = prepareForGetDP(P,node,elem,rnge_elec,rnge_gel,hdrInfo,elecNeeded,uniTag)
%
% Prepare to solve in getDP
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% October 2017

[dirname,baseFilename] = fileparts(P);
if isempty(dirname), dirname = pwd; end

% node = node + 0.5; already done right after mesh

indNode_elecElm = elem(find(elem(:,5) == 8),1:4);
X = zeros(size(indNode_elecElm,1),3);
for e = 1:size(indNode_elecElm,1), X(e,:) = mean ( node(indNode_elecElm(e,:),1:3) ); end
% figure; plot3(X(:,1),X(:,2),X(:,3),'r.');

label_elec = zeros(size(X,1),1);
for i = 1:size(X,1)
    offset = 0;
    while label_elec(i)==0
        for k = 1:length(rnge_elec)
            if ~isempty(rnge_elec{k}) && X(i,1)>rnge_elec{k}(2,1)-offset && X(i,1)<rnge_elec{k}(1,1)+offset && X(i,2)>rnge_elec{k}(2,2)-offset && X(i,2)<rnge_elec{k}(1,2)+offset && X(i,3)>rnge_elec{k}(2,3)-offset && X(i,3)<rnge_elec{k}(1,3)+offset
                label_elec(i) = k;
            end
        end
        offset = offset+0.1;
    end
    %     if mod(i,100) == 0
    %         fprintf('%f%% completed...\n',i*100/length(X))
    %     end
end
% figure; plot3(X(find(label_elec==1),1),X(find(label_elec==1),2),X(find(label_elec==1),3),'r.')

indNode_gelElm = elem(find(elem(:,5) == 7),1:4);
X = zeros(size(indNode_gelElm,1),3);
for e = 1:size(indNode_gelElm,1), X(e,:) = mean ( node(indNode_gelElm(e,:),1:3) ); end
% figure; plot3(X(:,1),X(:,2),X(:,3),'r.');

label_gel = zeros(size(X,1),1);
for i = 1:size(X,1)
    offset = 0;
    while label_gel(i)==0
        for k = 1:length(rnge_gel)
            if ~isempty(rnge_gel{k}) && X(i,1)>rnge_gel{k}(2,1)-offset && X(i,1)<rnge_gel{k}(1,1)+offset && X(i,2)>rnge_gel{k}(2,2)-offset && X(i,2)<rnge_gel{k}(1,2)+offset && X(i,3)>rnge_gel{k}(2,3)-offset && X(i,3)<rnge_gel{k}(1,3)+offset
                label_gel(i) = k;
            end
        end
        offset = offset+0.1;
    end
    %     if mod(i,100) == 0
    %         fprintf('%f%% completed...\n',i*100/length(X))
    %     end
end
% figure; plot3(X(find(label_gel==1),1),X(find(label_gel==1),2),X(find(label_gel==1),3),'r.')
save([dirname filesep baseFilename '_' uniTag '_elecMeshLabels.mat'],'label_elec','label_gel');

element_elecNeeded = cell(length(elecNeeded),1);
area_elecNeeded = zeros(length(elecNeeded),1);

resolution = mean(hdrInfo.pixdim);
% mean() here to handle anisotropic resolution; ugly. Maybe just
% resample MRI to isotropic in the very beginning?

warning('off','MATLAB:TriRep:PtsNotInTriWarnId');
for i=1:length(element_elecNeeded)
    
    if isempty(indNode_elecElm(label_elec==i,:))
        error(['Electrode ' elecNeeded{i} ' goes out of image boundary. ROAST cannot proceed without a properly placed electrode. Please expand the input MRI by specifying the ''zeroPadding'' option.']);
    end
        
    [faces_elec,verts_elec] = freeBoundary(TriRep(indNode_elecElm(label_elec==i,:),node(:,1:3)));
    [faces_gel,verts_gel] = freeBoundary(TriRep(indNode_gelElm(label_gel==i,:),node(:,1:3)));
    [~,iE,iG] = intersect(verts_elec,verts_gel,'rows');
    tempTag = ismember(faces_elec,iE);
    % faces_overlap = faces_elec(sum(tempTag,2)==3,:);
    faces_elecOuter = faces_elec(~(sum(tempTag,2)==3),:);
    [~,Loc] = ismember(verts_elec,node(:,1:3),'rows');
    element_elecNeeded{i} = Loc(faces_elecOuter);
    % calculate the surface area
    a = (verts_elec(faces_elecOuter(:, 2),:) - verts_elec(faces_elecOuter(:, 1),:))*resolution;
    b = (verts_elec(faces_elecOuter(:, 3),:) - verts_elec(faces_elecOuter(:, 1),:))*resolution;
    c = cross(a, b, 2);
    area_elecNeeded(i) = sum(0.5*sqrt(sum(c.^2, 2)));
    
end
save([dirname filesep baseFilename '_' uniTag '_usedElecArea.mat'],'area_elecNeeded');

disp('setting up boundary conditions...');

fid_in = fopen([dirname filesep baseFilename '_' uniTag '.msh']);
fid_out = fopen([dirname filesep baseFilename '_' uniTag '_ready.msh'],'w');

numOfTissue = length(unique(elem(:,5)));
while ~feof(fid_in)
    s = fgetl(fid_in);
    
    if strcmp(s,'$Elements')
        fprintf(fid_out,'%s\n',s);
        s = fgetl(fid_in);
        numOfElem = str2num(s);
        fprintf(fid_out,'%s\n',num2str(numOfElem+size(cell2mat(element_elecNeeded),1)));
    elseif strcmp(s,'$EndElements')
        ii = 0;
        for j=1:length(element_elecNeeded)
            for i=1:size(element_elecNeeded{j},1)
                
                fprintf(fid_out,'%s \n',[num2str(numOfElem+i+ii) ' 2 2 ' num2str(numOfTissue+j) ' ' num2str(numOfTissue+j) ' ' num2str(element_elecNeeded{j}(i,1)) ' ' num2str(element_elecNeeded{j}(i,2)) ' ' num2str(element_elecNeeded{j}(i,3))]);
                
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