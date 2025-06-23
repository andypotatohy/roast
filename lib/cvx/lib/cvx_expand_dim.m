function zx = cvx_expand_dim( z, dim, nx )

%CVX_EXPAND_DIM   Expands an N-D array along a specified dimension.
%   CVX_EXPAND_DIM( X, DIM, NX ) stacks NX copies of the matrix X along the
%   dimension NX. It is equivalent to CAT( DIM, X, X, ..., X ), where X is
%   repeated NX times.
%
%   This is an internal CVX function, and as such no checking is performed to
%   insure that the arguments are valid.

zdims = cell( 1, max( ndims(z), dim ) );
[ zdims{:} ] = deal( ':' );
zdims{dim} = ones( 1, nx );
zx = z( zdims{:} );

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

