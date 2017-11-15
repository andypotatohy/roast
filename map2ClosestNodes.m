function [closestDist,ind] = map2ClosestNodes(inputNodes,goalNodes)
% [closestDist,ind] = map2ClosestNodes(inputNodes,goalNodes)
%
% Note this function will give "out of memory" error if the input and goal
% point clounds are too big (>50K)

N_inputNodes = size(inputNodes,1);
N_goalNodes = size(goalNodes,1);

N_dim_input = size(inputNodes,2);
N_dim_goal = size(goalNodes,2);
% now supports any number of dimensions

if N_dim_input ~= N_dim_goal
    error('You cannot map points that are in spaces of different dimensions!')
end

temp = reshape(inputNodes',1,N_dim_input,N_inputNodes);

temp1 = repmat(temp,[N_goalNodes 1 1]);

temp2 = repmat(goalNodes,[1 1 N_inputNodes]);

dist = sqrt(sum((temp1 - temp2).^2,2));

dist = squeeze(dist);

if N_inputNodes == N_goalNodes && all(inputNodes(:) == goalNodes(:))
    dist(find(dist==0)) = inf;
end % if input and goal nodes are the same set of nodes

[closestDist,ind] = min(dist);  % can output the first argument, as a record of the minimal distance,
% so that can get an average distance after sweeping through all the
% cortical locations, and can be used as a performance measure for the
% Euclidean distance mapping