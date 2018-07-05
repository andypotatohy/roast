function opt = writeRoastLog(subject,opt)
% opt = writeRoastLog(subject,opt)
%
% Write simulation log.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

[dirname,baseFilename] = fileparts(subject);
if isempty(dirname), dirname = pwd; end

fid = fopen([dirname filesep baseFilename '_log'],'a');

if ~isempty(opt.uniqueTag)
    uniqueTag = opt.uniqueTag;
else
    %     uniqueTag = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss')); % for Matlab 2014b and later
    uniqueTag = char(datestr(now,30)); % for Matlab 2014a and earlier
    opt.uniqueTag = uniqueTag;
end

if ~exist([dirname filesep baseFilename '_' uniqueTag '_options.mat'],'file')
    save([dirname filesep baseFilename '_' uniqueTag '_options.mat'],'opt');
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

fprintf(fid,'zeroPadding:\t');
if opt.zeroPad>0
    fprintf(fid,'%d',opt.zeroPad);
else
    fprintf(fid,'none');
end

fprintf(fid,'\n\n');
fclose(fid);