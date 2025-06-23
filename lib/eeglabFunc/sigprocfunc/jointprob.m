% jointprob() - rejection of odd columns of a data array  using 
%              joint probability of the values in that column (and
%              using the probability distribution of all columns).
%
% Usage:
%   >>  [jp rej] = jointprob( signal );
%   >>  [jp rej] = jointprob( signal, threshold, jp, normalize, discret);
%  
%
% Inputs:
%   signal     - one dimensional column vector of data values, two 
%                dimensional column vector of values of size 
%                sweeps x frames or three dimensional array of size 
%                component x sweeps x frames. If three dimensional, 
%                all components are treated independently. 
%   threshold  - Absolute threshold. If normalization is used then the 
%                threshold is expressed in standard deviation of the
%                mean. 0 means no threshold.
%   jp         - pre-computed joint probability (only perform thresholding). 
%                Default is the empty array [].
%   normalize  - 0 = do not not normalize entropy. 1 = normalize entropy.
%                2 is 20% trimming (10% low and 10% high) proba. before 
%                normalizing. Default is 0.
%   discret    - discretization variable for calculation of the 
%                discrete probability density. Default is 1000 points. 
% 
% Outputs:
%   jp         - normalized joint probability  of the single trials 
%                (size component x sweeps)
%   rej        - rejected matrix (0 and 1, size comp x sweeps)
%
% Remark:
%   The exact values of joint-probability depend on the size of a time 
%   step and thus cannot be considered as absolute.
%
% See also: realproba()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [jp, rej] = jointprob( signal, threshold, oldjp, normalize, discret );

if nargin < 1
	help jointprob;
	return;
end;	
if nargin < 2
	threshold = 0;
end;	
if nargin < 3
	oldjp = [];
end;	
if nargin < 4
	normalize = 0;
end;	
if nargin < 5
	discret = 1000;
end;	

if size(signal,2) == 1 % transpose if necessary
	signal = signal';
end;

[nbchan pnts sweeps] = size(signal);
jp  = zeros(nbchan,sweeps);

if exist('oldjp') & ~isempty( oldjp ) % speed up the computation
	jp = oldjp;
else
	for rc = 1:nbchan

		% COMPUTE THE DENSITY FUNCTION
		% ----------------------------
		[ dataProba sortbox ] = realproba( signal(rc, :), discret );

		% compute all entropy
		% -------------------
		for index=1:sweeps
			datatmp = dataProba((index-1)*pnts+1:index*pnts);
			jp(rc, index) = - sum( log( datatmp ) ); 
			     % - sum( datatmp .* log( datatmp ) ); would be the entropy
		end;
	end;

	% normalize the last dimension
	% ----------------------------	
	if normalize
        tmpjp = jp;
        if normalize == 2,
            tmpjp = sort(jp);
            tmpjp = tmpjp(round(length(tmpjp)*0.1):end-round(length(tmpjp)*0.1));
        end;
        try, 
            switch ndims( signal )
             case 2,	jp = (jp-mean(tmpjp)) / std(tmpjp);
             case 3,	jp = (jp-mean(tmpjp,2)*ones(1,size(jp,2)))./ ...
                  (std(tmpjp,0,2)*ones(1,size(jp,2)));
            end;
        catch, error('Error while normalizing'); end;
	end;
end	

% reject
% ------	
if threshold ~= 0 
    if length(threshold) > 1
    	rej = (threshold(1) > jp) | (jp > threshold(2));
    else
    	rej = abs(jp) > threshold;
    end;
else
	rej = zeros(size(jp));
end;	

return;
