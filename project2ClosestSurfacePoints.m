function [distMetric_sorted,indexOnSurf] = project2ClosestSurfacePoints(points,surf,surfCenter)

vecP = points - repmat(surfCenter,size(points,1),1);
% vector connecting center to the point that is to be projected onto a surface
normVecP = repmat(sqrt(sum(vecP.^2,2)),1,size(vecP,2));
vecP = vecP./normVecP;
% for i=1:size(vecP,1), vecP(i,:) = vecP(i,:)/norm(vecP(i,:)); end
vecP = single(vecP);

vecS = surf - repmat(surfCenter,size(surf,1),1);
% vectors connecting center to each point on the surface
normVecS = repmat(sqrt(sum(vecS.^2,2)),1,size(vecS,2));
vecS = vecS./normVecS;
% for i=1:size(vecS,1), vecS(i,:) = vecS(i,:)/norm(vecS(i,:)); end
vecS = single(vecS);

temp = reshape(vecP',1,size(vecP,2),size(vecP,1));
vecP = repmat(temp,[size(vecS,1) 1 1]);
vecS = repmat(vecS,[1 1 size(vecP,3)]);
distMetric = squeeze(dot(vecP,vecS,2));
[distMetric_sorted,indexOnSurf] = sort(distMetric,'descend');

% distMetric_sorted = zeros(size(vecS,1),size(vecP,1));
% indexOnSurf = zeros(size(vecS,1),size(vecP,1));
% 
% for i=1:size(vecP,1)
%     distMetric = dot(repmat(vecP(i,:),size(vecS,1),1),vecS,2);
%     [distMetric_sorted(:,i),indexOnSurf(:,i)] = sort(distMetric,'descend');
% end