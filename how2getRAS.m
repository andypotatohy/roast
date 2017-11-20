function [permIn,permOut,flipTag] = how2getRAS(landmarks)

% Landmarks follow the order of: nasion, inion, right, left, front neck,
% and back neck.

% disp('adjusting the orientation of the head into RAS orientation...')
e1 = landmarks(3,:) - landmarks(4,:); % right-left
e1 = e1/norm(e1);

e2 = landmarks(1,:) - landmarks(2,:); % nasion-inion
e2 = e2/norm(e2);

e3 = landmarks(1,:) - landmarks(5,:); % nasion-front_neck
e3 = e3/norm(e3);
% detect the orientation of the head based on the anatomical landmarks

[~,perm1] = max(abs(e1)); [~,perm2] = max(abs(e2)); [~,perm3] = max(abs(e3));
flipTag = [sign(e1(perm1)) sign(e2(perm2)) sign(e3(perm3))];
% detect if the head is flipped or not in each direction compared to RAS system

permIn = [perm1,perm2,perm3]; % permutation order into RAS
[~,permOut] = sort(permIn); % inverse permutation order back to original orientation of the head
