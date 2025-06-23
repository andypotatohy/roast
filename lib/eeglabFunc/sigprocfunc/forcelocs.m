% forcelocs() - rotate location in 3-D so specified electrodes
%               match specified locations. 
%               CAUTION: Only for use on electrodes in
%               and remaining in the upper spherical hemisphere,
%               otherwise it will work improperly. Written primarily for
%               adjusting all electrodes homogenously with Cz.
%
% Usage:
%   >> chanlocs = forcelocs( chanlocs ); % pop-up window mode
%   >> chanlocs = forcelocs( chanlocs, loc1, loc2, ... );
% Example:
%   >> chanlocs = forcelocs( chanlocs, { 0.78, 'x', 'A1' }, { 0.023, 'x', ...
%   'B1','B2','Cz' } );
%
% Inputs:
%   chanlocs  - EEGLAB channel structure. See help readlocs()
%
% Optional inputs:
%   loc1      - cell array: { location, axis, channame1, channame2, .. } 
%               'location' is new cartesian coordinate of channame1 along 'axis'
%               'axis' is either
%                 'X'   New x-coordinate of mean of channame1, channame2,
%                       etc. Used to calculate the X-Z plane angle by
%                       which to rotate all channels.
%                       Note that all rotations are to the corresponding positive 
%                       Z-value, since theta=atan(z/x).
%                 'Y'   New x-coordinate of mean of channame1, channame2,
%                       etc.
%                 
%               'channame#'  Name of channel(s) to be rotated, as they appear in 
%                       chanlocs.label
%   loc2      - same as loc1
%
% Outputs:
%   chanlocs  - updated EEGLAB channel structure.
%
%
% Author: Arnaud Delorme, CNL / Salk Institute, 15 April 2003
%
% See also: readlocs()

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function [chanlocs,options] = forcelocs( chanlocs, varargin)
    
    NENTRY = 1; % number of lines in GUI
    FIELDS = { 'X' 'Y' };
    
    options = [];
    if nargin < 1
        help forcelocs;
        return;
    end;
    if nargin < 2
        geom = { [0.4 1 1 0.3] };
        uilist = { { 'style' 'text' 'string' 'X/Y value' 'tag' 'valstr' } ...
                   { 'style' 'text' 'string' 'Coordinate' } ...
                   { 'style' 'text' 'string' 'Electrode list' } ...
                   { } };
        for index = 1:NENTRY
            tag = [ 'c' int2str(index) ];
            geom = { geom{:}  [0.3 1 1 0.3] };
            uilist = { uilist{:},{ 'style' 'edit' 'string' fastif(index==1, '0','') }, ...
                       { 'style' 'listbox' 'string' 'X (rotate X-Z plane)|Y (rotate Y-Z plane)' ...
                         'callback' [ 'if get(gco, ''value'') == 1,' ...
                                      '     set(findobj(gcbf, ''tag'', ''valstr''), ''string'', ''Y value'');' ...
                                      'else set(findobj(gcbf, ''tag'', ''valstr''), ''string'', ''X value'');' ...
                                      'end;' ] }, ...
                       { 'style' 'edit' 'string'  fastif(index==1, 'Cz','') 'tag' tag }, ...
                       { 'style' 'pushbutton' 'string' 'Pick' ...
                         'callback', [ 'tmp3 = get(gcbf, ''userdata'');' ...
                                       '[tmp1 tmp2] = pop_chansel({tmp3.labels}, ''selectionmode'', ''single'');' ...
                                       'if ~isempty(tmp1) set(findobj(gcbf, ''tag'', ''' tag '''), ''string'', tmp2); end;' ...
                                       'clear tmp1 tmp2;' ] } };
        end;
        
        results = inputgui( geom, uilist, 'pophelp(''forcelocs'');', 'Force electrode location -- forcelocs()', chanlocs );
        if length(results) == 0, return; end;
        
        options = {};
        for index = 1:NENTRY
            tmpi = 3*(index-1)+1;
            if ~isempty(results{tmpi})
                tmpchans = parsetxt(results{tmpi+2});
                options = { options{:},{ str2num(results{tmpi}) FIELDS{results{tmpi+1}} tmpchans{:} }};
            end;
        end;    
    else 
        options = varargin;
    end;

    % scan all locations
    % ------------------
    channelnames = lower(strvcat({chanlocs.labels}));
    for index = 1:length(options)
        
        val   = options{index}{1};
        type  = options{index}{2};
        chans = getchans(options{index}(3:end), channelnames);

        % rotate X-Z plane 
        % ----------------
        if strcmpi(type, 'x')
            curx   = mean([ chanlocs(chans).X ]);
            curz   = mean([ chanlocs(chans).Z ]);
            newx = val;
            rotangle = solvesystem(curx, curz, newx);
            
            for chanind = 1:length(chanlocs)
                [chanlocs(chanind).X chanlocs(chanind).Z]= rotation(chanlocs(chanind).X, chanlocs(chanind).Z, rotangle);
            end;
            chanlocs = convertlocs(chanlocs, 'cart2all');
        end;
        
        % rotate Y-Z plane 
        % ----------------
        if strcmpi(type, 'y')
            cury   = mean([ chanlocs(chans).Y ]);
            curz   = mean([ chanlocs(chans).Z ]);
            newy = val;
            rotangle = solvesystem(cury, curz, newy);
            
            for chanind = 1:length(chanlocs)
                [chanlocs(chanind).Y chanlocs(chanind).Z]= rotation(chanlocs(chanind).Y, chanlocs(chanind).Z, rotangle);
            end;
            chanlocs = convertlocs(chanlocs, 'cart2all');
        end;
    
    end;
        

% get channel indices
% -------------------
function chanlist = getchans(chanliststr, channelnames);
    chanlist = [];
    for index = 1:length(chanliststr)
        i = strmatch (lower(chanliststr{index}), channelnames, 'exact');
        chanlist  = [chanlist i];
    end;

% function rotate coordinates
% ---------------------------
function [X,Y] = rotation(x,y,rotangle)
    X = real((x+j*y)*exp(j*rotangle));
    Y = imag((x+j*y)*exp(j*rotangle));
    
% function solvesyst
% ------------------
function theta = solvesystem(x,y,nx)
    % Original Solution
    %eq(1,:) = [x -y]; res(1) = nx;
    %eq(2,:) = [y x];  res(2) = sqrt(x^2+y^2-nx^2);
    %sol = eq\res';
    %theta = atan2(sol(2), sol(1));
    
    % simplified solution
    ny = sqrt(x^2+y^2-nx^2);
    ang1 = angle(x+j*y);
    ang2 = angle(nx+j*ny);
    theta = ang2-ang1;
    
    % Even simpler solution     Toby 03/05/2007
    % theta = atan(y/x);
    
