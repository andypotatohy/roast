% sbplot() - create axes in arbitrary subplot grid positions and sizes
%
% Usage:  >> axis_handle = sbplot(v,h,index)  
%         >> axis_handle = sbplot(v,h,[index1 index2])
%         >> axis_handle = sbplot(v,h,[index1 index2],axprop,..)
%         >> axis_handle = sbplot(v,h,[index1 index2],'ax',handle,axprop,..)
%
% Inputs:
%   v,h    - Integers giving the vertical and horizontal ranks of the tiling.
%   index  - Either a single subplot index, in which case the command 
%            is equivalent to subplot, or a two-element vector giving 
%            the indices of two corners of the sbplot() area according 
%            to subplot() convention (e.g., left-to-right, top-to-bottom).
%   axprop - Any axes property(s), e.g., >> sbplot(3,3,3,'color','w')
%   handle - Following keyword 'ax', sbplot tiles the given axes handle
%            instead of the whole figure
%
% Output:
%  axis_handle - matlab axis handle
%
% Note:
%   sbplot is essentially the same as the subplot command except that 
%   sbplot axes may span multiple tiles. Also, sbplot() will not erase 
%   underlying axes. 
%  
%  Examples: >> sbplot(3,3,6);plot(rand(1,10),'g');
%            >> sbplot(3,3,[7 2]);plot(rand(1,10),'r');
%            >> sbplot(8,7,47);plot(rand(1,10),'b'); 
%
% Authors: Colin Humphries, Arnaud Delorme & Scott Makeig, SCCN/INC/UCSD, La Jolla, June, 1998 

% Copyright (C) June 1998, Colin Humphries & Scott Makeig, SCCN/INC/UCSD, 
% scott@sccn.ucsd.edu
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

% reformatted by Scott Makeig, 6/10/98
% 12/22/00 test nargin<3 -sm
% 01/21/01 added (recursive) axes option 'ax' -sm
% 01-25-02 reformated help & licence -ad 

function [out] = sbplot(m,n,gridpos,varargin) % varargin is std. matlab arg list

if nargin<3
   error('            requires >=3 arguments');
end

if nargin>3 & strcmp(varargin{1},'ax')
   pos = get(varargin{2},'Position'); %  sbplot(3,1,[2 3]) -> 0.4111 0.1100 0.4939 0.815
   varargin = {varargin{3:end}};
else
   pos = get(gcf,'DefaultAxesPosition'); %  [0.1300    0.1100    0.7750    0.815]
end

Xpad = pos(1);      % lower-left distance from left side of figure
Ypad = pos(2);      % lower-left distance from bottom of figure

Xlen = pos(3);      % axes width
Ylen = pos(4);      % axes height

if n == 2
  xspace = Xlen*0.27/(n-0.27);      % xspace between axes as per subplot
else
  xspace = (0.9*Xlen)*0.27/(n-0.9*0.27); 
end

if m == 2
  yspace = Ylen*0.27/(m-0.27);      % yspace between axes as per subplot
else                                % WHY Xlen (.775) instead of Ylen (.815) ??
  yspace = (0.9*Ylen)*0.27/(m-0.9*0.27); 
end
  
xlength = (Xlen-xspace*(n-1))/n;  % axes width
ylength = (Ylen-yspace*(m-1))/m;  % axes height

% Convert tile indices to grid positions
if length(gridpos) == 1
  xgridpos(1) = mod(gridpos,n);   % grid position
  if xgridpos(1) == 0
    xgridpos(1) = n;                 
  end
  xgridpos(2) = 1;                   % grid length
  ygridpos(1) = m-ceil(gridpos/n)+1; % grid position
  ygridpos(2) = 1;                   % grid length
else
  xgridpos(1) = mod(gridpos(1),n);
  if xgridpos(1) == 0
    xgridpos(1) = n;
  end
  tmp = mod(gridpos(2),n);
  if tmp == 0
    tmp = n;
  end
  if tmp > xgridpos(1)
    xgridpos(2) = tmp-xgridpos(1)+1;
  else
    xgridpos(2) = xgridpos(1)-tmp+1;
    xgridpos(1) = tmp;
  end
  
  ygridpos(1) = m-ceil(gridpos(1)/n)+1;
  tmp = m-ceil(gridpos(2)/n)+1;
  if tmp > ygridpos(1)
    ygridpos(2) = tmp-ygridpos(1)+1;
  else 
    ygridpos(2) = ygridpos(1)-tmp+1;
    ygridpos(1) = tmp;
  end
end

% Calculate axes coordinates
position(1) = Xpad+xspace*(xgridpos(1)-1)+xlength*(xgridpos(1)-1);
position(2) = Ypad+yspace*(ygridpos(1)-1)+ylength*(ygridpos(1)-1)-0.03;
position(3) = xspace*(xgridpos(2)-1)+xlength*xgridpos(2);
position(4) = yspace*(ygridpos(2)-1)+ylength*ygridpos(2);

% Create new axes
ax = axes('Position',position,varargin{:});

% Output axes handle
if nargout > 0
  out = ax;         
end
