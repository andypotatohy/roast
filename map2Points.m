function [dist,indexOnGoalPoints]= map2Points(inputPoints,goalPoints,criterion,numOfPts)
% [dist,indexOnGoalPoints]= map2Points(inputPoints,goalPoints,criterion,numOfPts)
%
% Map from one point cloud to another point cloud, based on the Euclidean
% distance.
% 
% Note this function will give "out of memory" error if the input and goal
% point clounds are too big (>50K)

inputPoints = single(inputPoints);
goalPoints = single(goalPoints);

N_inputPoints = size(inputPoints,1);
N_goalPoints = size(goalPoints,1);

N_dim_input = size(inputPoints,2);
N_dim_goal = size(goalPoints,2);
% now supports any number of dimensions

if N_dim_input ~= N_dim_goal
    error('You cannot map points that are in spaces of different dimensions!')
end

temp = reshape(inputPoints',1,N_dim_input,N_inputPoints);

temp1 = repmat(temp,[N_goalPoints 1 1]);

temp2 = repmat(goalPoints,[1 1 N_inputPoints]);

dist = sqrt(sum((temp1 - temp2).^2,2));

dist = squeeze(dist);

[dist_sorted,ind_sortedDist] = sort(dist);

switch criterion
    case 'closest'
        dist = dist_sorted(1,:);
        indexOnGoalPoints = ind_sortedDist(1,:);
    case 'farthest'
        dist = dist_sorted(end,:);
        indexOnGoalPoints = ind_sortedDist(end,:);
    case 'closer'
        if isempty(numOfPts), error('You want to map to the first XX closest points, please specify XX in the 4th argument.'); end
        if numOfPts>N_goalPoints, error('Number of points exceed size of goal point cloud.'); end
        dist = dist_sorted(1:numOfPts,:);
        indexOnGoalPoints = ind_sortedDist(1:numOfPts,:);
    case 'farther'
        if isempty(numOfPts), error('You want to map to the first XX farthest points, please specify XX in the 4th argument.'); end
        if numOfPts>N_goalPoints, error('Number of points exceed size of goal point cloud.'); end
        dist = dist_sorted(end-numOfPts+1:end,:);
        indexOnGoalPoints = ind_sortedDist(end-numOfPts+1:end,:);
    otherwise
        error('Please specify either closest, farthest, closer or farther as the mapping criterion.');
end