function res = optimize_anon(p,t,A)
% res = optimize_anon(p,t,A)
%
% Objective function for the searching for the optimal orientation of the
% electric field at the target. Refer to the following paper for details:
% J. P. Dmochowski, A. Datta, Y. Huang, J. D. Richardson, M. Bikson, J. Fridriksson,
% L. C. Parra, Targeted transcranial direct current stimulation for rehabilitation
% after stroke, NeuroImage, Vol. 75, July 2013, pp. 12–19
%
% Variable to be optimized:
% t, an n-by-2 matrix, where n is the number of target ROIs. Each row in t corresponds
% to each target ROI, representing the electric field orientation by a UNIT vector
% in spherical coordinates, which is made up of:
% t(n,1) -- azimuth; t(n,2) -- elevation.
%
% Other parts of this function are almost the same as the gateway function
% for optimization (optimize.m). Refer to the following
% paper for details: J.P. Dmochowski, A. Datta, M. Bikson, Y. Su, L.C. Parra
% Optimized multi-electrode stimulation increases focality and intensity at target
% J. Neural Eng., 8 (2011), p. 046011
%
% This implement accepts any number of targeting ROIs.
%
% optimize_prepare.m should be called before calling this function.
%
% Jacek P. Dmochowski, 2012
% Yu (Andy) Huang, October 2014
% Yu (Andy) Huang, January 2017

if size(t,2)~=2
    error('Orientation variable not in correct format. Please provide the orientation in n-by-2 matrix, where n is the number of target ROIs, and 2 columns representing azimuth and elevation.')
end

% A = p.A;
numOfTargets = p.numOfTargets;
Nlocs = p.Nlocs;
% node_distances = p.node_distances;
% sorted_nodes = p.sorted_nodes;
target_nodes = p.target_nodes;

optType = p.optType;
elecNum = p.elecNum;
I_max = p.I_max;

w = p.w;
U = p.U; S = p.S; V = p.V;

Ed = zeros(numOfTargets,3);
for n = 1:numOfTargets
    Ed(n,:) = [cos(t(n,1))*sin(t(n,2)) sin(t(n,1))*sin(t(n,2)) cos(t(n,2))]; % spherical to cartesian coordinates
end

a = p.a;
xd = zeros(3*Nlocs,1);
for n = 1:numOfTargets
    xd(target_nodes{n}) = a*Ed(n,1);
    xd(target_nodes{n}+Nlocs) = a*Ed(n,2);
    xd(target_nodes{n}+2*Nlocs) = a*Ed(n,3);
end

% CORE ALGORITHM
[xopt,~,status] = optimize_currents(A,xd,I_max,w,target_nodes,optType,U,S,V,elecNum,0);
if strcmp(status,'Failed')
    warning('Warn:convert',...
        'Inner optimization FAILED!!\n Program will continue but results may be INACCURATE!\n');
    % else
    %     fprintf('Inner optimization COMPLETED successfully!\n')
end

% Output of this objective function
if ~isempty(strfind(optType,'wls'))
    res = norm(sqrt(w).*(xd - xopt));
% targetint = zeros(numOfTargets,1);
%     for n=1:numOfTargets
%         targetint(n) = dot( Ed(n,:) , mean ( [ xopt(target_nodes{n}) , xopt(target_nodes{n}+Nlocs) , xopt(target_nodes{n}+2*Nlocs) ], 1 ) );
%     end
%     res = -sum(targetint);
elseif ~isempty(strfind(optType,'lcmv'))
    res = norm(xopt);
elseif ~isempty(strfind(optType,'max-l1'))
    targetint = zeros(numOfTargets,1);
    for n=1:numOfTargets
        targetint(n) = dot( Ed(n,:) , mean ( [ xopt(target_nodes{n}) , xopt(target_nodes{n}+Nlocs) , xopt(target_nodes{n}+2*Nlocs) ], 1 ) );
    end
    res = -sum(targetint);
end