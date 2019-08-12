function [x,cvx_status] = max_l1per(f,ub,solElecNum,verbose)
% Accepts any number of targeting ROIs.
% solElecNum is the desired number of anodes+cathodes in the optimal
% solution/montage.
% ANDY 2014-10-27
% ANDY 2017-01-30

n = size(f,1);

if verbose
    cvx_begin
else
    cvx_begin quiet
end
              variable x(n);
              maximize( sum(f'*x) );
               subject to
                 norm(x,1)+abs(sum(x)) <= ub;
                 norm([x; sum(x)],inf) <= ub/solElecNum;
                 %          sum(x)==0;
    cvx_end