function isNewOpt = isNewOptions(optNew,optOld)

isNewOpt = 0;

if ~strcmp(optNew.configTxt,optOld.configTxt)
    isNewOpt = 1;
    return
end

capNew = [any(strcmp(optNew.elecPara(1).capType,{'1020','1010','1005'})) any(strcmp(optNew.elecPara(1).capType,{'biosemi','Biosemi','bioSemi','BioSemi','BIOSEMI'})) strcmp(optNew.elecPara(1).capType,'none')];
capOld = [any(strcmp(optOld.elecPara(1).capType,{'1020','1010','1005'})) any(strcmp(optOld.elecPara(1).capType,{'biosemi','Biosemi','bioSemi','BioSemi','BIOSEMI'})) strcmp(optOld.elecPara(1).capType,'none')];
if any(capNew~=capOld)
    isNewOpt = 1;
    return
end

if length(optNew.elecPara) > length(optOld.elecPara)
    temp = repmat(optOld.elecPara,length(optNew.elecPara),1);
    if size(optOld.elecPara.elecSize,1)>1
        for i=1:length(temp), temp(i).elecSize = temp(i).elecSize(i,:); end
    end
    if size(optOld.elecPara.elecOri,1)>1
        for i=1:length(temp), temp(i).elecOri = temp(i).elecOri(i,:); end
    end
    optOld.elecPara = temp;
elseif length(optNew.elecPara) < length(optOld.elecPara)
    temp = repmat(optNew.elecPara,length(optOld.elecPara),1);
    if size(optNew.elecPara.elecSize,1)>1
        for i=1:length(temp), temp(i).elecSize = temp(i).elecSize(i,:); end
    end
    if size(optNew.elecPara.elecOri,1)>1
        for i=1:length(temp), temp(i).elecOri = temp(i).elecOri(i,:); end
    end
    optNew.elecPara = temp;
end

for i=1:length(optNew.elecPara)
    
    elecNew = [any(strcmp(optNew.elecPara(i).elecType,{'disc','Disc','DISC'})) ...
        any(strcmp(optNew.elecPara(i).elecType,{'pad','Pad','PAD'})) ...
        any(strcmp(optNew.elecPara(i).elecType,{'ring','Ring','RING'}))];
    elecOld = [any(strcmp(optOld.elecPara(i).elecType,{'disc','Disc','DISC'})) ...
        any(strcmp(optOld.elecPara(i).elecType,{'pad','Pad','PAD'})) ...
        any(strcmp(optOld.elecPara(i).elecType,{'ring','Ring','RING'}))];
    if any(elecNew~=elecOld), isNewOpt = 1; return; end
    
    if size(optNew.elecPara(i).elecSize,1)>1 && size(optOld.elecPara(i).elecSize,1)==1
        temp = optOld.elecPara(i).elecSize;
        optOld.elecPara(i).elecSize = repmat(temp,size(optNew.elecPara(i).elecSize,1),1);
    elseif size(optNew.elecPara(i).elecSize,1)==1 && size(optOld.elecPara(i).elecSize,1)>1
        temp = optNew.elecPara(i).elecSize;
        optNew.elecPara(i).elecSize = repmat(temp,size(optOld.elecPara(i).elecSize,1),1);
    end
    for j=1:size(optNew.elecPara(i).elecSize,1)
        if any(optNew.elecPara(i).elecSize(j,:) ~= optOld.elecPara(i).elecSize(j,:))
           isNewOpt = 1;
           return
        end
    end
    
    if find(elecNew)==2
        if ischar(optNew.elecPara(i).elecOri) && ischar(optOld.elecPara(i).elecOri)
            if ~strcmp(optNew.elecPara(i).elecOri,optOld.elecPara(i).elecOri)
                isNewOpt = 1;
                return
            end
        elseif ~ischar(optNew.elecPara(i).elecOri) && ~ischar(optOld.elecPara(i).elecOri)
            if size(optNew.elecPara(i).elecOri,1)>1 && size(optOld.elecPara(i).elecOri,1)==1
                temp = optOld.elecPara(i).elecOri;
                optOld.elecPara(i).elecOri = repmat(temp,size(optNew.elecPara(i).elecOri,1),1);
            elseif size(optNew.elecPara(i).elecOri,1)==1 && size(optOld.elecPara(i).elecOri,1)>1
                temp = optNew.elecPara(i).elecOri;
                optNew.elecPara(i).elecOri = repmat(temp,size(optOld.elecPara(i).elecOri,1),1);
            end
            for j=1:size(optNew.elecPara(i).elecOri,1)
                if any(optNew.elecPara(i).elecOri(j,:) ~= optOld.elecPara(i).elecOri(j,:))
                   isNewOpt = 1;
                   return
                end
            end
        else
            isNewOpt = 1; return
        end
    end
    
end

if ~(isempty(optNew.T2) && isempty(optOld.T2))
    if ~strcmp(optNew.T2,optOld.T2)
        isNewOpt = 1;
        return
    end
end

if optNew.meshOpt.radbound~=optOld.meshOpt.radbound
    isNewOpt = 1;
    return
end

if optNew.meshOpt.angbound~=optOld.meshOpt.angbound
    isNewOpt = 1;
    return
end

if optNew.meshOpt.distbound~=optOld.meshOpt.distbound
    isNewOpt = 1;
    return
end

if optNew.meshOpt.reratio~=optOld.meshOpt.reratio
    isNewOpt = 1;
    return
end

if optNew.meshOpt.maxvol~=optOld.meshOpt.maxvol
    isNewOpt = 1;
    return
end