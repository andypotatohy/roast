function p = optimize_prepare(p,A,locs)
% p = optimize_prepare(p,A,locs)
%
% Preprocessing function for running the optimization.
%
% This function calculates target info and algorithm parameters (except
% field orientation at targets). Singular value decomposition of matrix A
% will be done here to save time.
%
% After this function, the gateway function for optimization (optimize.m)
% should be called; or if one wants to search for optimal orientation at
% targets, optimize_anon.m should be called, and optimize.m will be called
% afterwards.
%
% Refer to the following papers for details:
% J.P. Dmochowski, A. Datta, M. Bikson, Y. Su, L.C. Parra
% Optimized multi-electrode stimulation increases focality and intensity at target
% J. Neural Eng., 8 (2011), p. 046011
% and
% J. P. Dmochowski, A. Datta, Y. Huang, J. D. Richardson, M. Bikson, J. Fridriksson,
% L. C. Parra, Targeted transcranial direct current stimulation for rehabilitation
% after stroke, NeuroImage, Vol. 75, July 2013, pp. 12–19
%
% This implement accepts any number of targeting ROIs.
%
% Yu (Andy) Huang, October 2014
% Yu (Andy) Huang, January 2017

% A = p.A;
% locs = p.locs;
Nnodes = p.Nnodes;

targetCoord = p.targetCoord;
numOfTargets = p.numOfTargets;

% node_distances = zeros(Nnodes,numOfTargets);
% sorted_nodes = zeros(Nnodes,numOfTargets);
distances_to_target = zeros(Nnodes,numOfTargets);
for n = 1:numOfTargets
    tmp = locs - repmat(targetCoord(n,:),size(locs,1),1);
    distances_to_target(:,n) = sqrt(sum(tmp.*tmp,2));
    %     [node_distances(:,n),sorted_nodes(:,n)] = sort(distances_to_target(:,n));
end

targetRadius = p.targetRadius;
target_nodes = cell(numOfTargets,1);
for n = 1:numOfTargets
    target_nodes{n} = find(distances_to_target(:,n)<targetRadius);
    if isempty(target_nodes{n})
        error('No nodes found near target. Please increase the value of ''targetRadius''.');
    end
end

p.target_nodes = target_nodes;
target_nodes = cell2mat(target_nodes);
nontarget_nodes = setdiff(1:Nnodes,target_nodes);

optType = p.optType;

fprintf('Preparing to optimize...this may take a moment...\n');

if ~isempty(strfind(optType,'wls'))
    k = p.k;
    Ntarget = length(target_nodes);
    Nnontarget = Nnodes-Ntarget;
    wp = Nnodes/Ntarget*(k/(k+1));
    wn = Ntarget/Nnontarget*(wp/k);
else
    wp = 1;
    wn = 0;
end

w = zeros(Nnodes,1);
w(target_nodes) = wp;
w(nontarget_nodes) = wn;
w = repmat(w,3,1);
p.w = w;

if strcmp(optType,'wls-l1') || strcmp(optType,'wls-l1per') || strcmp(optType,'wls-l1penalty')
    [U,S,V] = svd(repmat(sqrt(w),1,size(A,2)).*A, 0);
elseif strcmp(optType,'lcmv-l1') || strcmp(optType,'lcmv-l1per')
    [U,S,V] = svd(A,0);
else
    U = []; S = []; V = [];
end
p.U = U; p.S = S; p.V = V;

% fprintf('Preprocessing done! Ready for optimization!\n');