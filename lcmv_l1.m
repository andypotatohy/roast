function [x,cvx_status] = lcmv_l1(A,C,f,ub,verbose)
% Accepts any number of targeting ROIs.
% ANDY 2014-10-27
% ANDY 2017-01-30

n = size(A,2);

if verbose
    cvx_begin
else
    cvx_begin quiet
end
              variable x(n);
              minimize( norm(A*x,2) );
               subject to
                  norm([x;-sum(x)],1) <= 2*ub;
                  C*x==f;
    cvx_end