function y = cvx_id( x )
y = apply( @cvx_id, x );
switch class( y ),
    case 'struct',
        y = struct2cell( y );
        y = max( [ -Inf, y{:} ] );
    case 'cell',
        y = max( [ -Inf, y{:} ] );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.


