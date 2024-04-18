function renameSPMres(src,tar)
% renameSPMres(src,tar)
%
% (c) Andrew Birnbaum, Parra Lab at CCNY
%     Yu (Andy) Huang
% April 2024

[dir,srcName] = fileparts(src);
if isempty(dir), dir = pwd; end
[~,tarName] = fileparts(tar);

for t = 1:6
    movefile([dir filesep 'c' num2str(t) srcName '.nii'], ...
        [dir filesep 'c' num2str(t) tarName '.nii']);
end

movefile([dir filesep srcName '_rmask.mat'], ...
    [dir filesep tarName '_rmask.mat']);

movefile([dir filesep srcName '_seg8.mat'], ...
    [dir filesep tarName '_seg8.mat']);
