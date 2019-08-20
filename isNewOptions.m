function isNewOpt = isNewOptions(optNew,optOld,fun)
% isNewOpt = isNewOptions(optNew,optOld,fun)
%
% Compare two options to see if they are different.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
% August 2019 callable by roast_target()

isNewOpt = 0;

switch fun
    
    case 'roast'
        
        if ~strcmp(optNew.configTxt,optOld.configTxt)
            isNewOpt = 1;
            return
        end
        
        capNew = [any(strcmp(optNew.elecPara(1).capType,{'1020','1010','1005'})) strcmpi(optNew.elecPara(1).capType,'biosemi') strcmp(optNew.elecPara(1).capType,'none')];
        capOld = [any(strcmp(optOld.elecPara(1).capType,{'1020','1010','1005'})) strcmpi(optOld.elecPara(1).capType,'biosemi') strcmp(optOld.elecPara(1).capType,'none')];
        if any(capNew~=capOld)
            isNewOpt = 1;
            return
        end
        
        % below complicated logics for checking if new options for electrode sizes
        % and orientations (for pad), because 'elecSize' and 'elecOri' can be very
        % flexible
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
            
            elecNew = [strcmpi(optNew.elecPara(i).elecType,'disc') ...
                strcmpi(optNew.elecPara(i).elecType,'pad') ...
                strcmpi(optNew.elecPara(i).elecType,'ring')];
            elecOld = [strcmpi(optOld.elecPara(i).elecType,'disc') ...
                strcmpi(optOld.elecPara(i).elecType,'pad') ...
                strcmpi(optOld.elecPara(i).elecType,'ring')];
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
                    if ~strcmpi(optNew.elecPara(i).elecOri,optOld.elecPara(i).elecOri)
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
                else % needs improvement here
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
        
        if optNew.conductivities.white~=optOld.conductivities.white
            isNewOpt = 1;
            return
        end
        
        if optNew.conductivities.gray~=optOld.conductivities.gray
            isNewOpt = 1;
            return
        end
        
        if optNew.conductivities.csf~=optOld.conductivities.csf
            isNewOpt = 1;
            return
        end
        
        if optNew.conductivities.bone~=optOld.conductivities.bone
            isNewOpt = 1;
            return
        end
        
        if optNew.conductivities.skin~=optOld.conductivities.skin
            isNewOpt = 1;
            return
        end
        
        if optNew.conductivities.air~=optOld.conductivities.air
            isNewOpt = 1;
            return
        end
        
        if length(optNew.conductivities.gel)~=length(optOld.conductivities.gel)
            isNewOpt = 1;
            return
        else
            if any(optNew.conductivities.gel~=optOld.conductivities.gel)
                isNewOpt = 1;
                return
            end
        end
        
        if length(optNew.conductivities.electrode)~=length(optOld.conductivities.electrode)
            isNewOpt = 1;
            return
        else
            if any(optNew.conductivities.electrode~=optOld.conductivities.electrode)
                isNewOpt = 1;
                return
            end
        end
        
        if optNew.resamp ~= optOld.resamp
            isNewOpt = 1;
            return
        end
        
        if optNew.zeroPad ~= optOld.zeroPad
            isNewOpt = 1;
            return
        end
        
    case 'target'
        
        if ~strcmp(optNew.roastTag,optOld.roastTag)
            isNewOpt = 1;
            return
        end
        
        if size(optNew.targetCoord,1) ~= size(optOld.targetCoord,1)
            isNewOpt = 1;
            return
        end
        
        [isDone,indInOldRun] = ismember(optNew.targetCoord,optOld.targetCoord,'rows');
        if ~all(isDone), isNewOpt = 1; return; end
        % indInOldRun ensures target order will not affect determining if
        % it's a new run
        
        u0 = optOld.u0(indInOldRun,:);
        if any(optNew.u0(:) ~= u0(:))
            isNewOpt = 1;
            return
        end
        
        if ~iscell(optOld.orient)
            orient = optOld.orient(indInOldRun,:);
        else
            orient = optOld.orient(indInOldRun);
        end
        % below handles 'optimal' and 'radial-in', as these two share the
        % same u0
        if all(optNew.u0(:) == u0(:))
            if (iscell(optNew.orient) && strcmpi(optNew.orient{1},'optimal')) ...
                    && (~iscell(orient) || (iscell(orient) && ~strcmpi(orient{1},'optimal')))
                isNewOpt = 1;
                return
            end
            if (iscell(orient) && strcmpi(orient{1},'optimal')) ...
                    && (~iscell(optNew.orient) || (iscell(optNew.orient) && ~strcmpi(optNew.orient{1},'optimal')))
                isNewOpt = 1;
                return
            end
        end
        
        if ~strcmpi(optNew.optType,optOld.optType)
            isNewOpt = 1;
            return
        end
        
        if optNew.targetRadius ~= optOld.targetRadius
            isNewOpt = 1;
            return
        end
        
        if strcmpi(optNew.optType,'max-l1per')
            if optNew.elecNum ~= optOld.elecNum
                isNewOpt = 1;
                return
            end
        end
        
        if ~isempty(strfind(lower(optNew.optType),'wls'))
            if optNew.k ~= optOld.k
                isNewOpt = 1;
                return
            end
        end
        
        if strcmpi(optNew.optType,'lcmv-l1') || strcmpi(optNew.optType,'lcmv-l1per')
            if optNew.a ~= optOld.a
                isNewOpt = 1;
                return
            end
        end
        
end