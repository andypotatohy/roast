function [x_opt,s_opt,status] = optimize_currents(A,d,S_max,w,tar_nodes,method,U,S,V,solElecNum,verbose)
% [x_opt,s_opt,status] = optimize_currents(A,d,S_max,w,tar_nodes,method,U,S,V,solElecNum,verbose)
%
% Core implementation for the tDCS targeting technique. Refer to the following
% paper for details: J.P. Dmochowski, A. Datta, M. Bikson, Y. Su, L.C. Parra
% Optimized multi-electrode stimulation increases focality and intensity at target
% J. Neural Eng., 8 (2011), p. 046011
%
% INPUTS:
% A: transfer matrix (or lead field in EEG community, gain matrix in EIT community),
% a 3N-by-M matrix where N is the number of volume mesh nodes and M is the
% number of electrodes (excluding reference electrode)
%
% d: desired electric field (3N vector where N is the number of nodes)
%
% S_max: maximum allowable current density
%
% w: weights used for weighted least square method, also used for indexing
% target areas
%
% tar_nodes: indices of target areas, stored separately for each area/ROI,
% for calculating mean value for each ROI
%
% solElecNum: desired number of anodes+cathodes in the optimal solution (i.e.
% electrode montage). ONLY applies to max intensity with L1 constraint on
% individual electrode, i.e., max-l1per.
%
% verbose: with screen output (1, default) or without screen output (0)
%
% OUTPUTS:
% x_opt: optimized electric field (3*N vector where N is the number of nodes)
%
% s_opt: optimal applied current densities (M vector where M is the number of
% electrodes, excluding reference electrode)
%
% status: a status variable indicating if the optimization is successful or
% not
%
% This implement accepts any number of targeting ROIs.
%
% Jacek P. Dmochowski, 2011
% Yu (Andy) Huang, October 2014
% Yu (Andy) Huang, January 2017

if nargin < 11
    verbose = 1;
end

M = size(A,2);
Nlocs = size(A,1)/3;
numOfROI = length(tar_nodes);


if strcmp(method,'unconstrained-wls') % unconstrained weighted least squares
    sqrtw = sqrt(w);
    s_opt = (repmat(sqrtw,1,M).*A) \ (sqrtw.*d);
    x_opt = A*s_opt;
    status = 'Solved';
elseif strcmp(method,'wls-l1') % weighted least squares with L1 constraint
    sqrtw = sqrt(w);
    [s_opt,status] = ls_l1(S*V',U'*(sqrtw.*d),2*S_max,verbose);
    x_opt = A*s_opt;
elseif strcmp(method,'wls-l1per') % weighted least squares with L1 constraint on individual electrode
    sqrtw = sqrt(w);
    [s_opt,status] = ls_l1per(S*V',U'*(sqrtw.*d),S_max,verbose);
    x_opt = A*s_opt;
elseif strcmp(method,'wls-l1penalty') % weighted least squares with L1 constraint as a penalty in the cost function
                                      % individual constraint on injected current will be imposed when optimization is done
    sqrtw = sqrt(w);
    [s_opt,status] = ls_l1asPenalty(S*V',U'*(sqrtw.*d),0.01*2*S_max,verbose);
    % 0.01 is the best lambda so far (it should be relied on data)
    x_opt = A*s_opt;
    
elseif strcmp(method,'unconstrained-lcmv') % unconstrained LCMV
    C = zeros(3*numOfROI,M);
    f = zeros(3*numOfROI,1);
    for n = 1:numOfROI
        C(((n-1)*3+1):n*3,:) = [mean(A(tar_nodes{n},:),1);mean(A(tar_nodes{n}+Nlocs,:),1);mean(A(tar_nodes{n}+2*Nlocs,:),1);];
        f(((n-1)*3+1):n*3) = [mean(d(tar_nodes{n}));mean(d(tar_nodes{n}+Nlocs));mean(d(tar_nodes{n}+2*Nlocs));];
    end
    ATAi = inv(A'*A);
    s_opt = ATAi*C'*inv(C*ATAi*C')*f;
    x_opt = A*s_opt;
    status = 'Solved';
elseif strcmp(method,'lcmv-l1') % LCMV with L1 constraint
    C = zeros(3*numOfROI,M);
    f = zeros(3*numOfROI,1);
    for n = 1:numOfROI
        C(((n-1)*3+1):n*3,:) = [mean(A(tar_nodes{n},:),1);mean(A(tar_nodes{n}+Nlocs,:),1);mean(A(tar_nodes{n}+2*Nlocs,:),1);];
        f(((n-1)*3+1):n*3) = [mean(d(tar_nodes{n}));mean(d(tar_nodes{n}+Nlocs));mean(d(tar_nodes{n}+2*Nlocs));];
    end
    status = '';
    while strcmp(status,'Solved')~=1
        %         fprintf('CVX optimization infeasible under LCMV with L1 constraint...reducing target intensity...\n')
        [s_opt,status] = lcmv_l1(S*V',C,f,2*S_max,verbose);
        f = f/1.25;
    end
    if all(abs(f)<1e-5)
        s_opt = nan(M,1);
        %         error('CVX optimization infeasible for this target using LCMV with L1 constraint! Consider targeting a nearby location.')
    end
    x_opt = A*s_opt;
elseif strcmp(method,'lcmv-l1per') % LCMV with L1 constraint on individual electrode
    C = zeros(3*numOfROI,M);
    f = zeros(3*numOfROI,1);
    for n = 1:numOfROI
        C(((n-1)*3+1):n*3,:) = [mean(A(tar_nodes{n},:),1);mean(A(tar_nodes{n}+Nlocs,:),1);mean(A(tar_nodes{n}+2*Nlocs,:),1);];
        f(((n-1)*3+1):n*3) = [mean(d(tar_nodes{n}));mean(d(tar_nodes{n}+Nlocs));mean(d(tar_nodes{n}+2*Nlocs));];
    end
    status = '';
    while strcmp(status,'Solved')~=1
        %         fprintf('CVX optimization infeasible under LCMV with L1 constraint on individual electrode...reducing target intensity...\n')
        [s_opt,status] = lcmv_l1per(S*V',C,f,S_max,verbose);
        f = f/1.25;
    end
    if all(abs(f)<1e-5)
        s_opt = nan(M,1);
        %         error('CVX optimization infeasible for this target using LCMV with L1 constraint on individual electrode! Consider targeting a nearby location.')
    end
    x_opt = A*s_opt;
    
elseif strcmp(method,'max-l1') % maximum intensity in desired direction with L1 constraint
    
    %     indx = find(w);
    %     C = A(indx,:);
    %     f = d(indx);
    %     Cf = C'*f;    
    Cf = zeros(M,numOfROI);
    % C = zeros(3*numOfROI,M);
    for n = 1:numOfROI
        %         indx1 = ;
        %         indx2 = ;
        %         indx3 = ;
        %     Cx = A(indx1,:);
        %     Cy = A(indx2,:);
        %     Cz = A(indx3,:);
        Cx = mean(A(tar_nodes{n},:),1);
        Cy = mean(A(tar_nodes{n}+Nlocs,:),1);
        Cz = mean(A(tar_nodes{n}+2*Nlocs,:),1);
        C = [Cx;Cy;Cz];
        %     C(((n-1)*3+1):n*3,:) = [Cx;Cy;Cz];
        %     f = [d(indx1);d(indx2);d(indx3)];
        f = [mean(d(tar_nodes{n}));mean(d(tar_nodes{n}+Nlocs));mean(d(tar_nodes{n}+2*Nlocs))];
        Cf(:,n) = C'*f;
    end
    
    [s_opt,status] = max_l1(Cf,2*S_max,verbose);
    % [s_opt,status] = max_l1norm(Cf,2*S_max);
    
    x_opt = A*s_opt;
elseif strcmp(method,'max-l1per') % maximum intensity in desired direction with L1 constraint on individual electrode
    %     indx = find(w);
    %     C = A(indx,:);
    %     f = d(indx);
    %     Cf = C'*f;
    
    Cf = zeros(M,numOfROI);
    for n = 1:numOfROI
        Cx = mean(A(tar_nodes{n},:),1);
        Cy = mean(A(tar_nodes{n}+Nlocs,:),1);
        Cz = mean(A(tar_nodes{n}+2*Nlocs,:),1);
        C = [Cx;Cy;Cz];
        f = [mean(d(tar_nodes{n}));mean(d(tar_nodes{n}+Nlocs));mean(d(tar_nodes{n}+2*Nlocs))];
        Cf(:,n) = C'*f;
    end
    
    [s_opt,status] = max_l1per(Cf,2*S_max,solElecNum,verbose);
    x_opt = A*s_opt;
end