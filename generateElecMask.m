function vol_elec = generateElecMask(elec_allCoord,coordRange,elec,doWarn)
% vol_elec = generateElecMask(elec_allCoord,coordRange,elec,doWarn)
% 
% Convert point clouds of electrodes into 3D masks.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

vol_elec = zeros(coordRange);
% vol_elecLabel = zeros(coordRange);
% rnge = cell(length(elec_allCoord),1);

for i = 1:length(elec_allCoord)
    temp = elec_allCoord{i};
    if isempty(temp)
        error(['Electrode ' elec{i} ' goes out of image boundary. ROAST cannot proceed without a properly placed electrode. Please expand the input MRI by specifying the ''zeroPadding'' option.']);
    else
        ind = find(temp(:,1)>0 & temp(:,1)<=coordRange(1)...
            & temp(:,2)>0 & temp(:,2)<=coordRange(2)...
            & temp(:,3)>0 & temp(:,3)<=coordRange(3));
        if isempty(ind)
            error(['Electrode ' elec{i} ' goes out of image boundary. ROAST cannot proceed without a properly placed electrode. Please expand the input MRI by specifying the ''zeroPadding'' option.']);
%             elec_allCoord{i} = [];
%             rnge{i} = [];
        else
            if length(ind)<size(temp,1)
                if doWarn
                    warning(['Part of the electrode ' elec{i} ' goes out of image boundary. ROAST can continue but results may not be accurate. It is recommended that you expand the input MRI by specifying the ''zeroPadding'' option.']);
                end
            end
            temp = temp(ind,:);
%             can detect elec overlap here
%             vol_elec(sub2ind(size(vol_elec),temp(:,1),temp(:,2),temp(:,3))) = 1;
            vol_elec(sub2ind(size(vol_elec),temp(:,1),temp(:,2),temp(:,3))) = i; % individual mask/color for each elec & gel
%             vol_elecLabel(sub2ind(size(vol_elec),temp(:,1),temp(:,2),temp(:,3))) = i;
%             elec_allCoord{i} = temp;
%             rnge{i} = [max(temp);min(temp)];
        end
    end
end

% allCoord = cell2mat(elec_allCoord);
% vol_elec(sub2ind(size(vol_elec),allCoord(:,1),allCoord(:,2),allCoord(:,3)))=1;