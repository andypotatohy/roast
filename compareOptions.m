function isNewOpt = compareOptions(optNew,optOld)

isNewOpt = 0;

if ~strcmp(optNew.configTxt,optOld.configTxt)
    isNewOpt = 1;
    return
end

if length(optNew.elecPara) == length(optOld.elecPara)
    capNew = [any(strcmp(optNew.elecPara(1).capType,{'1020','1010','1005'})) any(strcmp(optNew.elecPara(1).capType,{'biosemi','Biosemi','bioSemi','BioSemi','BIOSEMI'}))];
    capOld = [any(strcmp(optOld.elecPara(1).capType,{'1020','1010','1005'})) any(strcmp(optOld.elecPara(1).capType,{'biosemi','Biosemi','bioSemi','BioSemi','BIOSEMI'}))];
    if capNew~=capOld
        isNewOpt = 1;
        return
    end
    for i=1:length(optNew.elecPara)
        optNew.elecPara(i).elecType
    end
else
   
    if length(optNew.elecPara) > length(optOld.elecPara)
        
    else
        
    end
    
end
% if length(optNew.elecPara)==1 && length(optOld.elecPara)==1
%     ~any(strcmp(capType,{'1020','1010','1005','biosemi','Biosemi','bioSemi','BioSemi','BIOSEMI'}))
%     
% elseif length(optNew.elecPara)>1 && length(optOld.elecPara)>1
%     
% elseif length(optNew.elecPara)>1 && length(optOld.elecPara)==1
%     
% elseif length(optNew.elecPara)==1 && length(optOld.elecPara)>1
%     
% end
 
% t2
% 
% meshopt