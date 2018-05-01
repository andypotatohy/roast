function [vol_elec,rnge,elec_allCoord] = generateElecMask(elec_allCoord,coordRange,elec,doWarn)
% [vol_elec,rnge,elec_allCoord] = generateElecMask(elec_allCoord,coordRange,elec,doWarn)
% 
% Convert point clouds of electrodes into 3D masks.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

vol_elec = zeros(coordRange);
rnge = cell(length(elec_allCoord),1);

for i = 1:length(elec_allCoord)
    temp = elec_allCoord{i};
    if ~isempty(temp)
        ind = find(temp(:,1)>0 & temp(:,1)<=coordRange(1)...
            & temp(:,2)>0 & temp(:,2)<=coordRange(2)...
            & temp(:,3)>0 & temp(:,3)<=coordRange(3));
        if isempty(ind)
            error(['Electrode ' elec{i} ' goes out of image boundary. ROAST cannot solve without a properly placed electrode. You can do zero-padding to the MRI (run zeroPadding()) before running ROAST.']);
            elec_allCoord{i} = [];
            rnge{i} = [];
        else
            if length(ind)<size(temp,1)
                if doWarn
                    warning(['Part of the electrode ' elec{i} ' goes out of image boundary. ROAST can continue but results may not be accurate. It is recommended that you do zero-padding to the MRI (run zeroPadding()) before running ROAST.']);
                end
            end
            temp = temp(ind,:);
            elec_allCoord{i} = temp;
            rnge{i} = [max(temp);min(temp)];
        end
    end
end

allCoord = cell2mat(elec_allCoord);
vol_elec(sub2ind(size(vol_elec),allCoord(:,1),allCoord(:,2),allCoord(:,3)))=1;