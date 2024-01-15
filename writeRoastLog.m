function opt = writeRoastLog(subject,opt,fun)
% opt = writeRoastLog(subject,opt,fun)
%
% Write logs for roast() and roast_target().
% The log file for roast() ends with '_roastLog' and the log file for
% roast_target() ends with '_targetLog'.
%
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018
% August 2019 callable by roast_target()

[dirname,baseFilename] = fileparts(subject);
if isempty(dirname), dirname = pwd; end

switch fun
    
    case 'roast'
        
        fid = fopen([dirname filesep baseFilename '_roastLog'],'a');
        
        if ~isempty(opt.uniqueTag)
            uniqueTag = opt.uniqueTag;
        else
            %     uniqueTag = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss')); % for Matlab 2014b and later
            uniqueTag = char(datestr(now,30)); % for Matlab 2014a and earlier
            opt.uniqueTag = uniqueTag;
        end
        
        if ~exist([dirname filesep baseFilename '_' uniqueTag '_roastOptions.mat'],'file')
            save([dirname filesep baseFilename '_' uniqueTag '_roastOptions.mat'],'opt');
        else
            error('You''re about to run a simulation using options that you never used before, but forgot to use a new tag for it. ROAST will get confused when managing different simulations with the same tag. Please use a new tag.');
        end
        
        fprintf(fid,'%s:\n',uniqueTag);
        fprintf(fid,'recipe:\t%s\n',opt.configTxt);
        fprintf(fid,'capType:\t%s\n',opt.elecPara(1).capType);
        
        fprintf(fid,'elecType:\t');
        for i=1:length(opt.elecPara)
            fprintf(fid,'%s\t',opt.elecPara(i).elecType);
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'elecSize:\t');
        if length(opt.elecPara)==1
            for i=1:size(opt.elecPara.elecSize,1)
                fprintf(fid,'[');
                for j=1:size(opt.elecPara.elecSize,2)
                    fprintf(fid,'%.1f',opt.elecPara.elecSize(i,j));
                    if j<size(opt.elecPara.elecSize,2), fprintf(fid,','); end
                end
                fprintf(fid,']');
                if i<size(opt.elecPara.elecSize,1), fprintf(fid,'; '); end
            end
        else
            for i=1:length(opt.elecPara)
                fprintf(fid,'[');
                for j=1:size(opt.elecPara(i).elecSize,2)
                    fprintf(fid,'%.1f',opt.elecPara(i).elecSize(j));
                    if j<size(opt.elecPara(i).elecSize,2), fprintf(fid,','); end
                end
                fprintf(fid,']');
                if i<length(opt.elecPara), fprintf(fid,'; '); end
            end
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'elecOri:\t');
        if length(opt.elecPara)==1
            if isempty(opt.elecPara.elecOri)
                fprintf(fid,'N/A');
            elseif ischar(opt.elecPara.elecOri)
                fprintf(fid,'%s',opt.elecPara.elecOri);
            else
                for i=1:size(opt.elecPara.elecOri,1)
                    fprintf(fid,'[');
                    for j=1:size(opt.elecPara.elecOri,2)
                        fprintf(fid,'%.4f',opt.elecPara.elecOri(i,j));
                        if j<size(opt.elecPara.elecOri,2), fprintf(fid,','); end
                    end
                    fprintf(fid,']');
                    if i<size(opt.elecPara.elecOri,1), fprintf(fid,'; '); end
                end
            end
        else
            for i=1:length(opt.elecPara)
                if isempty(opt.elecPara(i).elecOri)
                    fprintf(fid,'N/A');
                elseif ischar(opt.elecPara(i).elecOri)
                    fprintf(fid,'%s',opt.elecPara(i).elecOri);
                else
                    fprintf(fid,'[');
                    for j=1:size(opt.elecPara(i).elecOri,2)
                        fprintf(fid,'%.4f',opt.elecPara(i).elecOri(j));
                        if j<size(opt.elecPara(i).elecOri,2), fprintf(fid,','); end
                    end
                    fprintf(fid,']');
                end
                if i<length(opt.elecPara), fprintf(fid,'; '); end
            end
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'T2:\t');
        if ~isempty(opt.T2)
            fprintf(fid,'%s',opt.T2);
        else
            fprintf(fid,'none');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'meshOpt:\t');
        fprintf(fid,'radbound: %d; angbound: %d; distbound: %.1f; reratio: %d; maxvol: %d',...
            opt.meshOpt.radbound,opt.meshOpt.angbound,opt.meshOpt.distbound,opt.meshOpt.reratio,opt.meshOpt.maxvol);
        fprintf(fid,'\n');
        
        fprintf(fid,'conductivities:\t');
        fprintf(fid,'white: %.3f; gray: %.3f; CSF: %.3f; bone: %.3f; skin: %.3f; air: %.1e; ',...
            opt.conductivities.white,opt.conductivities.gray,opt.conductivities.csf,opt.conductivities.bone,opt.conductivities.skin,opt.conductivities.air);
        fprintf(fid,'gel: ');
        for i=1:length(opt.conductivities.gel), fprintf(fid,'%.3f; ',opt.conductivities.gel(i)); end
        fprintf(fid,'electrode: ');
        for i=1:length(opt.conductivities.electrode), fprintf(fid,'%.1e; ',opt.conductivities.electrode(i)); end
        fprintf(fid,'\n');
        
        fprintf(fid,'reSampling:\t');
        if opt.resamp
            fprintf(fid,'on');
        else
            fprintf(fid,'off');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'multipriors:\t');
        if opt.multipriors
            fprintf(fid,'on');
        else
            fprintf(fid,'off');
        end
        fprintf(fid,'\n');

        fprintf(fid,'zeroPadding:\t');
        if opt.zeroPad>0
            fprintf(fid,'%d',opt.zeroPad);
        else
            fprintf(fid,'none');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'original MRI in RAS?:\t');
        if opt.isNonRAS
            fprintf(fid,'no, and is re-oriented into RAS');
        else
            fprintf(fid,'yes');
        end
        
        fprintf(fid,'\n\n');
        fclose(fid);
        
    case 'target'
        
        fid = fopen([dirname filesep baseFilename '_targetLog'],'a');
        
        if ~isempty(opt.uniqueTag)
            uniqueTag = opt.uniqueTag;
        else
            %     uniqueTag = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss')); % for Matlab 2014b and later
            uniqueTag = char(datestr(now,30)); % for Matlab 2014a and earlier
            opt.uniqueTag = uniqueTag;
        end
        
        if ~exist([dirname filesep baseFilename '_' uniqueTag '_targetOptions.mat'],'file')
            save([dirname filesep baseFilename '_' uniqueTag '_targetOptions.mat'],'opt');
        else
            error('You''re about to run a targeting using options that you never used before, but forgot to use a new tag for it. ROAST-TARGET will get confused when managing different targetings with the same tag. Please use a new tag.');
        end
        
        fprintf(fid,'%s:\n',uniqueTag);
        fprintf(fid,'ROAST simulation tag:\t%s\n',opt.roastTag);
        
        fprintf(fid,'targetCoord (in MNI space):\t');
        if ~isempty(opt.targetCoordMNI)
            for i=1:size(opt.targetCoordMNI,1)
                fprintf(fid,'[');
                for j=1:size(opt.targetCoordMNI,2)
                    fprintf(fid,'%.1f',opt.targetCoordMNI(i,j));
                    if j<size(opt.targetCoordMNI,2), fprintf(fid,','); end
                end
                fprintf(fid,']');
                if i<size(opt.targetCoordMNI,1), fprintf(fid,'; '); end
            end
        else
            fprintf(fid,'not provided');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'targetCoord (in original MRI voxel space):\t');
        if ~isempty(opt.targetCoordOriginal)
            for i=1:size(opt.targetCoordOriginal,1)
                fprintf(fid,'[');
                for j=1:size(opt.targetCoordOriginal,2)
                    fprintf(fid,'%d',opt.targetCoordOriginal(i,j));
                    if j<size(opt.targetCoordOriginal,2), fprintf(fid,','); end
                end
                fprintf(fid,']');
                if i<size(opt.targetCoordOriginal,1), fprintf(fid,'; '); end
            end
        else
            fprintf(fid,'not provided');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'targetCoord (in model voxel space):\t');
        for i=1:size(opt.targetCoord,1)
            fprintf(fid,'[');
            for j=1:size(opt.targetCoord,2)
                fprintf(fid,'%d',opt.targetCoord(i,j));
                if j<size(opt.targetCoord,2), fprintf(fid,','); end
            end
            fprintf(fid,']');
            if i<size(opt.targetCoord,1), fprintf(fid,'; '); end
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'desired intensity at each target (in V/m, only for max-focality):\t');
        if ~isempty(strfind(lower(opt.optType),'wls')) || ~isempty(strfind(lower(opt.optType),'lcmv'))
            fprintf(fid,'%.3f',opt.desiredIntensity);
        else
            fprintf(fid,'N/A');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'desired orientation at each target:\t');
        if ~iscell(opt.orient)
            for i=1:size(opt.orient,1)
                fprintf(fid,'[');
                for j=1:size(opt.orient,2)
                    fprintf(fid,'%.1f',opt.orient(i,j));
                    if j<size(opt.orient,2), fprintf(fid,','); end
                end
                fprintf(fid,']');
                if i<size(opt.orient,1), fprintf(fid,'; '); end
            end
        else
            for i=1:length(opt.orient)
                if ischar(opt.orient{i})
                    fprintf(fid,'%s',opt.orient{i});
                else
                    fprintf(fid,'[');
                    for j=1:size(opt.orient{i},2)
                        fprintf(fid,'%.4f',opt.orient{i}(j));
                        if j<size(opt.orient{i},2), fprintf(fid,','); end
                    end
                    fprintf(fid,']');
                end
                if i<length(opt.orient), fprintf(fid,'; '); end
            end
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'targeting algorithm:\t%s\n',opt.optType);
        
        fprintf(fid,'number of electrodes (only for max-l1per):\t');
        if ~isempty(opt.elecNum)
            fprintf(fid,'%d',opt.elecNum);
        else
            fprintf(fid,'N/A');
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'target radius (in mm):\t%d\n',opt.targetRadius);
        
        fprintf(fid,'k (only for wls):\t');
        if ~isempty(opt.k)
            fprintf(fid,'%.3f',opt.k);
        else
            fprintf(fid,'N/A');
        end
        fprintf(fid,'\n');
        
        fclose(fid);
        
    case 'target-results'
        
        fid = fopen([dirname filesep baseFilename '_targetLog'],'a');
        
        fprintf(fid,'optimal montage:\t%s\n',opt.montageTxt);
        
        fprintf(fid,'E-field magnitude at each target (in V/m):\t');
        for i=1:length(opt.targetMag)
            fprintf(fid,'%.2f',opt.targetMag(i));
            if i<length(opt.targetMag), fprintf(fid,','); end
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'E-field intensity along the desired orientation at each target (in V/m):\t');
        for i=1:length(opt.targetInt)
            fprintf(fid,'%.2f',opt.targetInt(i));
            if i<length(opt.targetInt), fprintf(fid,','); end
        end
        fprintf(fid,'\n');
        
        fprintf(fid,'E-field focality at each target (in cm):\t');
        for i=1:length(opt.targetMagFoc)
            fprintf(fid,'%.2f',opt.targetMagFoc(i));
            if i<length(opt.targetMagFoc), fprintf(fid,','); end
        end
        
        fprintf(fid,'\n\n');
        fclose(fid);
        
end