function seg8Name(P,P2,multipriors)

[dir, baseFile, ~] = fileparts(P);
[dir2, baseFile2, ~] = fileparts(P2);
if ~multipriors
    for t = 1:6
        copyfile([dir filesep 'c' num2str(t) baseFile '.nii'], ...
            [dir2 filesep 'c' num2str(t) baseFile2 '.nii']);
    end
    
    copyfile([dir filesep baseFile '_rmask.mat'], ...
        [dir2 filesep baseFile2 '_rmask.mat']);
    
    copyfile([dir filesep baseFile '_seg8.mat'], ...
        [dir2 filesep baseFile2 '_seg8.mat']);
else 
    copyfile([dir filesep baseFile '_seg8.mat'], ...
    [dir2 filesep baseFile2 '_seg8.mat']);
end
