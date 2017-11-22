function [neck_coord,neck_center]= placeNeckElec(scalp,scalp_surface,landmarks)
%
% Landmarks follow the order of: nasion, inion, right, left, front neck,
% and back neck.

% nasion = landmarks(1,:);
% inion = landmarks(2,:);
% right = landmarks(3,:);
% left = landmarks(4,:);
front_neck = landmarks(5,:);
back_neck = landmarks(6,:);

neck_center = (front_neck+back_neck)/2;

neck_elec = [front_neck;
    back_neck;
    neck_center(1)-round(size(scalp,1)/2) neck_center(2) neck_center(3);
    neck_center(1)+round(size(scalp,1)/2) neck_center(2) neck_center(3)];

idx = zeros(size(neck_elec,1),1);
[cosineAngle,indOnScalpSurf] = project2ClosestSurfacePoints(neck_elec,scalp_surface,neck_center);
for i = 1:length(idx)
%     testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > max(cosineAngle(:,i))*0.99993,i),:);
    testPts = scalp_surface(indOnScalpSurf(cosineAngle(:,i) > prctile(cosineAngle(:,i),99.99),i),:);
    [~,indFarthestOnTestPts] = map2Points(neck_center,testPts,'farthest');
    idx(i) = indOnScalpSurf(indFarthestOnTestPts,i);
    % Find the only point on the outer surface of the scalp for each electrode,
    % i.e., the exact coordinates for each electrode on the scalp surface
end

neck_coord = scalp_surface(idx,:);