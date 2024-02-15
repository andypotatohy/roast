function seg8Name(P, P2)

[dir, baseFile, ~] = fileparts(P);
[dir2, baseFile2, ~] = fileparts(P2);

for t = 1:6
    movefile([dir filesep 'c' num2str(t) baseFile '.nii'], ...
        [dir2 filesep 'c' num2str(t) baseFile2 '.nii']);
end

movefile([dir filesep baseFile '_rmask.mat'], ...
    [dir2 filesep baseFile2 '_rmask.mat']);

movefile([dir filesep baseFile '_seg8.mat'], ...
    [dir2 filesep baseFile2 '_seg8.mat']);
