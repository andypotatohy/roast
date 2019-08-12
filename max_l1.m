function [x,cvx_status] = max_l1(f,ub,verbose)
% Accepts any number of targeting ROIs.
% ANDY 2014-10-27
% ANDY 2017-01-30

n = size(f,1);

if verbose
    cvx_begin
else
    cvx_begin quiet
end
              variable x(n);
%               maximize( f'*x );
              maximize( sum(f'*x) );
               subject to
                 norm(x,1)+abs(sum(x)) <= ub;
%                  abs(sum(x)) <= ub;
%                  abs(x)<=ub;
% norm(x,1)+abs(sum(x)) <= ub;
% norm([x; sum(x)],inf) <= ub/4;
                 %          sum(x)==0;
    cvx_end