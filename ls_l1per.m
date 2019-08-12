function [x,cvx_status] = ls_l1per(A,d,ub,verbose)
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
              minimize( norm(A*x-d,2) );
               subject to
                 abs(sum(x))<= ub;
                 abs(x)<=ub;
    cvx_end