% writelocs() - write a file containing channel location, type and gain information
%             
% Usage:
%   >> writelocs( chanstruct, filename );
%   >> writelocs( chanstruct, filename, 'key', 'val' );
%
% Inputs:
%   chanstruct - EEG.chanlocs data structure returned by readlocs() containing
%                channel location, type and gain information.
%   filename   - File name for saving channel location, type and gain information
%
% Optional inputs:
%   'filetype'  - ['loc'|'sph'|'sfp'|'xyz'|'polhemus'|'besa'|'chanedit'|'custom'] 
%                 Type of the file to write. By default the file type is indicated 
%                 by the file extension. 
%                  'loc' - An EEGLAB 2-D polar coordinates channel locations file 
%                          Coordinates are theta and radius (see definitions below).
%                  'sph' - A Matlab spherical coordinates file (Note: spherical
%                          coordinates used by Matlab functions are different 
%                          from spherical coordinates used in BESA - see below).
%                  'sfp' - EGI cartesian coordinates (not Matlab cartesian - see below).
%                  'xyz' - MATLAB/EEGLAB cartesian coordinates (Not EGI cartesian; 
%                          z is toward nose; y is toward left ear; z is toward vertex).
%                  'polhemus' or 'polhemusx' - Polhemus electrode location file recorded with 
%                          'X' on sensor pointing to subject (see below and readelp()).
%                  'polhemusy' - Polhemus electrode location file recorded with 
%                          'Y' on sensor pointing to subject (see below and readelp()).
%                  'besa' - BESA'(.elp') spherical coordinate file. (Not MATLAB spherical
%                           - see below).
%                  'chanedit' - EEGLAB channel location files created by pop_chanedit().
%                  'custom' - Ascii files with columns in user-defined 'format' (see below).
%   'format'    - [cell array] Format of a 'custom' channel location file (see above).
%                          Default if no file type is defined. The cell array contains
%                          labels defining the meaning of each column of the input file.
%                           'channum'   [positive integer] channel number 
%                           'labels'    [string] channel name (no spaces)
%                           'theta'     [real degrees] 2-D angle in polar coordinates. 
%                                       positive => rotating from nose (0) toward left ear 
%                           'radius'    [real] radius in 2-D polar coords (0.5 is disk limits)
%                           'X'         [real] Matlab-cartesian X coordinate (to nose)
%                           'Y'         [real] Matlab-cartesian Y coordinate (to left ear)
%                           'Z'         [real] Matlab-cartesian Z coordinate (to vertex)
%                           '-X','-Y','-Z' Matlab-cartesian coordinates pointing away from above
%                           'sph_theta' [real degrees] Matlab spherical horizontal angle. 
%                                       positive => rotating from nose (0) toward left ear.
%                           'sph_phi'   [real degrees] Matlab spherical elevation angle;
%                                       positive => rotating from horizontal (0) upwards.
%                           'sph_radius' [real] distance from head center (unused) 
%                           'sph_phi_besa' [real degrees] BESA phi angle from vertical. 
%                                       positive => rotating from vertex (0) towards right ear.
%                           'sph_theta_besa' [real degrees] BESA theta horiz/azimuthal angle. 
%                                       positive => rotating from right ear (0) toward nose.
%     The input file may also contain other channel information fields
%                           'type'      channel type: 'EEG', 'MEG', 'EMG', 'ECG', others ...
%                           'calib'     [real near 1.0] channel calibration value.
%                           'gain'      [real > 1] channel gain. 
%                           'custom1'   custom field #1.
%                           'custom2', 'custom3', 'custom4' more custom fields.
%   'unicoord'     - ['on'|'off'] Uniformize all coordinates. Default 'on'.
%   'header'       - ['on'|'off'] Add a header comment line with the name of each column.
%                             Comment lines begin with '%'. Default is 'off'.
%   'customheader' - [string] Add a custom header at the beginning of the file and
%                             preceded by '%'.  If used with 'header' set to 'on', 
%                             the column names will be insterted after the custom header.
%   'elecind' - [integer array] Indices of channels to export. 
%                             Default is all channels.
%
% Note: for file formats, see readlocs() help  (>> help readlocs)
%
% Author: Arnaud Delorme, Salk Institute, 16 Dec 2002
%
% See also: readlocs(), readelp()

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 28 Feb 2002
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

function writelocs( chans, filename, varargin ); 

if nargin < 2
	help writelocs;
	return;
end;

% get infos from readlocs
% -----------------------
%[listtype formatinfo listcolformat formatskip] = readlocs('getinfos');
[chanformat listcolformat] = readlocs('getinfos');
indformat  = [];
for index = 1:length(chanformat), 
    if ~isstr(chanformat(index).importformat)
        indformat = [ indformat index ];
    end;
    if isempty(chanformat(index).skipline), chanformat(index).skipline = 0; end;
end;
listtype   = { chanformat(indformat).type };
formatinfo = { chanformat(indformat).importformat };
formatskip = [ chanformat(indformat).skipline ];

g = finputcheck( varargin, ...
                 { 'filetype'	  'string'	 listtype 			'loc';
                   'header'       'string'   { 'on','off' } 	'off';
                   'customheader' 'string'   [] 					'';
                   'elecind'      'integer'  [1 Inf]				[];
                   'unicoord'     'string'   { 'on','off' } 	'on'; 
                   'format'		  'cell'	 []					{} }, 'writelocs');
if isstr(g), error(g); end;  

if strcmpi(g.unicoord, 'on')
    disp('Uniformizing coordinates');
    chans = convertlocs(chans, 'auto', 'verbose', 'off');
end;

% select channels
% ---------------
if ~isempty(g.elecind)
	chans = chans(g.elecind);
end;

% finding types of input
% ----------------------
if isempty(g.format)
   indexformat = strmatch(lower(g.filetype), listtype, 'exact');
   g.format = formatinfo{indexformat};
   g.skipline = formatskip(indexformat);
else 
   g.skipline = 0;   
end;

% creating file
% -------------
fid = fopen(filename, 'w');

% exporting header
% ----------------
if ~isempty(g.customheader)
    allstrs = cellstr(g.customheader);
    for index=1:length(allstrs)
        fprintf(fid, '%s\n', allstrs{index});
    end;
end;
if  strcmpi(g.header, 'on') | g.skipline == 2
   for index=1:length(g.format)
      fprintf(fid, '%8s\t', g.format{index});
   end;
   fprintf(fid, '\n');
   for index=1:length(g.format)
      fprintf(fid, '%8s\t', char(ones(1,8)*45));
   end;
   fprintf(fid, '\n');
end;
if g.skipline == 1
   fprintf(fid, '%d\n', length(chans));
end;         

% writing infos
% -------------
for indexchan = 1:length(chans)
   for index=1:length(g.format)
      [str, mult] = checkformat(g.format{index});
      if strcmpi(str, 'channum')
         fprintf(fid, '%d', indexchan);
      else
         if ~isfield(chans, str)
            error([ 'Non-existant field: ''' str '''' ]);
         end;
         eval( [ 'chanval = chans(indexchan).' str ';' ] );
         if   isstr(chanval), fprintf(fid, '%8s', chanval);
         else   	
             if abs(mult*chanval) > 1E-10
                 fprintf(fid, '%8s', num2str(mult*chanval,5));
             else
                 fprintf(fid, '%8s', '0');
             end;
         end;
      end;
      if index ~= length(g.format)
         fprintf(fid, '\t');
      end;         
   end;
   fprintf(fid, '\n');
end;
fclose(fid);

return;

% check field format
% ------------------
function [str, mult] = checkformat(str)
	mult = 1;
	if strcmpi(str, 'labels'), str = lower(str); return; end;
	if strcmpi(str, 'channum'), str = lower(str); return; end;
	if strcmpi(str, 'theta'), str = lower(str); return; end;
	if strcmpi(str, 'radius'), str = lower(str); return; end;
	if strcmpi(str, 'sph_theta'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi'), str = lower(str); return; end;
	if strcmpi(str, 'sph_radius'), str = lower(str); return; end;
	if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi_besa'), str = lower(str); return; end;
	if strcmpi(str, 'gain'), str = lower(str); return; end;
	if strcmpi(str, 'calib'), str = lower(str); return; end;
	if strcmpi(str, 'type') , str = lower(str); return; end;
	if strcmpi(str, 'X'), str = upper(str); return; end;
	if strcmpi(str, 'Y'), str = upper(str); return; end;
	if strcmpi(str, 'Z'), str = upper(str); return; end;
	if strcmpi(str, '-X'), str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Y'), str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Z'), str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, 'custum1'), return; end;
	if strcmpi(str, 'custum2'), return; end;
	if strcmpi(str, 'custum3'), return; end;
	if strcmpi(str, 'custum4'), return; end;
   error(['writelocs: undefined field ''' str '''']);
   
