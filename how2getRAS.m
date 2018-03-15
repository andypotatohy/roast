function [permIn,permOut,flipTagInner,flipTagOutter] = how2getRAS(landmarks)

% Landmarks follow the order of: nasion, inion, right, left, front neck,
% and back neck.

nasion = landmarks(1,:);
inion = landmarks(2,:);
right = landmarks(3,:);
left = landmarks(4,:);
front_neck = landmarks(5,:);
% back_neck = landmarks(6,:);

% disp('adjusting the orientation of the head into RAS orientation...')
e1 = right - left;
e1 = e1/norm(e1);

e2 = nasion - inion;
e2 = e2/norm(e2);

e3 = nasion - front_neck;
e3 = e3/norm(e3);
% detect the orientation of the head based on the anatomical landmarks

[~,perm1] = max(abs(e1)); [~,perm2] = max(abs(e2)); [~,perm3] = max(abs(e3));
flipTagInner = [sign(e1(perm1)) sign(e2(perm2)) sign(e3(perm3))];
% detect if the head is flipped or not in each direction compared to RAS system

permIn = [perm1,perm2,perm3]; % permutation order into RAS
[~,permOut] = sort(permIn); % inverse permutation order back to original orientation of the head
flipTagOutter = flipTagInner(permOut);