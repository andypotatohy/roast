function r = optimize(p,A)
% r = optimize(p,A)
%
% Gateway function for the tDCS targeting technique. Refer to the following
% paper for details: J.P. Dmochowski, A. Datta, M. Bikson, Y. Su, L.C. Parra
% Optimized multi-electrode stimulation increases focality and intensity at target
% J. Neural Eng., 8 (2011), p. 046011
%
% This implement accepts any number of targeting ROIs.
%
% optimize_prepare.m should be called before calling this function.
%
% Jacek P. Dmochowski, 2011
% Yu (Andy) Huang, October 2014
% Yu (Andy) Huang, January 2017

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

Ed = p.u;
a = p.a;
xd = zeros(3*Nlocs,1);
for n = 1:numOfTargets
    xd(target_nodes{n}) = a*Ed(n,1);
    xd(target_nodes{n}+Nlocs) = a*Ed(n,2);
    xd(target_nodes{n}+2*Nlocs) = a*Ed(n,3);
end

% CORE ALGORITHM
fprintf('============================\nPerforming optimization...\n============================\n')
[xopt,sopt,status] = optimize_currents(A,xd,I_max,w,target_nodes,optType,U,S,V,elecNum,0);
if strcmp(status,'Failed')
    warning('Warn:convert',...
        'Optimization FAILED!!\n Program will continue but results may be INACCURATE!\n');
    % else
    %     fprintf('\n\nOptimization COMPLETED successfully!\n\n')
end

% OUTPUT RESULTS
% xoptmag = sqrt(xopt(1:Nlocs).^2+xopt(Nlocs+1:2*Nlocs).^2+xopt(2*Nlocs+1:3*Nlocs).^2);

% directivity = zeros(Nlocs,numOfTargets);
% crad = zeros(numOfTargets,1);
% for n = 1:numOfTargets
%     directivity(:,n) = cumsum( xoptmag(sorted_nodes(:,n)) ) / sum(xoptmag);
%     tmp = find(directivity(:,n)>0.5);
%     if ~isempty(tmp)
%         crad(n) = node_distances(tmp(1),n);
%     else
%         crad(n) = inf;
%     end
% end

r.sopt = sopt;
% r.xopt = xopt;
% r.xoptmag = xoptmag;

% r.directivity = directivity;
% r.crad = crad;

% r.targetintraw = xoptmag(sorted_nodes(1,:));
% r.targetint = dot(Ed,reshape([xopt(sorted_nodes(1,:)); xopt(sorted_nodes(1,:)+Nlocs);xopt(sorted_nodes(1,:)+2*Nlocs)],numOfTargets,3),2); % intensity in specified direction

% if exist('target_nodes','var')
targetintraw = zeros(numOfTargets,1);
targetint = zeros(numOfTargets,1);
for n=1:numOfTargets
    targetintraw(n) =  norm ( mean( [ xopt(target_nodes{n}) , xopt(target_nodes{n}+Nlocs) , xopt(target_nodes{n}+2*Nlocs) ], 1 ) );
    targetint(n) = dot( Ed(n,:) , mean ( [ xopt(target_nodes{n}) , xopt(target_nodes{n}+Nlocs) , xopt(target_nodes{n}+2*Nlocs) ], 1 ) );
end
r.targetintraw = targetintraw;
r.targetint = targetint;
% end