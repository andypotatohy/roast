function y = pow_pos( x, p )

%POW_POS   Power of positive part.
%   POW_POS(X,P) = POS(X).^P = MAX(X,0).^P.
%   Both P and X must be real, and P must be greater than or equal to 1.
%
%   Disciplined convex programming information:
%       POW_POS(X,P) is convex and nondecreasing in X; so when used in CVX
%       expressions, X must be convex. P must be constant, real, and
%       greater than or equal to 1.

error( nargchk( 2, 2, nargin ) );
if ~isnumeric( x ) || ~isreal( x ) || ~isnumeric( p ) || ~isreal( p ),
    error( 'Arguments must be real.' );
elseif any( p(:) <= 1 ),
    error( 'Second argument must be greater than or equal to 1.\nFor other exponents, use POW_P instead.', 1 ); %#ok
end
y = max(x,0).^p;

% Copyright 2012 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
