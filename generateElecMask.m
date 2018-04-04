function [vol_elec,rnge,elec_allCoord] = generateElecMask(elec_allCoord,coordRange,elec)

vol_elec = zeros(coordRange);
rnge = cell(length(elec_allCoord),1);

for i = 1:length(elec_allCoord)
    temp = elec_allCoord{i};
    if ~isempty(temp)
        ind = find(temp(:,1)>0 & temp(:,1)<=coordRange(1)...
            & temp(:,2)>0 & temp(:,2)<=coordRange(2)...
            & temp(:,3)>0 & temp(:,3)<=coordRange(3));
        if isempty(ind)
            warning(['Electrode ' elec{i} ' goes out of image boundary. you can do zero-padding to the MRI before running roast.']);
            elec_allCoord{i} = [];
            rnge{i} = [];
        else
            if length(ind)<size(temp,1)
                warning(['Part of the electrode ' elec{i} ' goes out of image boundary. you can do zero-padding to the MRI before running roast.']);
            end
            temp = temp(ind,:);
            elec_allCoord{i} = temp;
            rnge{i} = [max(temp);min(temp)];
        end
    end
end

allCoord = cell2mat(elec_allCoord);
vol_elec(sub2ind(size(vol_elec),allCoord(:,1),allCoord(:,2),allCoord(:,3)))=1;