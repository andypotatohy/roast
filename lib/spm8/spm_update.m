function spm_update(update)
% Check (and install) SPM8 updates from the FIL server
% FORMAT spm_update
% This function will connect itself to the FIL server, compare the
% version number of the updates with the one of the SPM installation 
% currently in the MATLAB path and will display the outcome.
%
% FORMAT spm_update(update)
% Invoking this function with any input parameter will do the same as
% above but will also attempt to download and install the updates.
%__________________________________________________________________________
% Copyright (C) 2010-2015 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_update.m 6312 2015-01-23 13:01:57Z guillaume $

url = 'http://www.fil.ion.ucl.ac.uk/spm/download/spm8_updates/';

if ~nargin
    update = false;
else
    update = true;
end

[s,sts] = urlread(url);
if ~sts, error('Cannot access the SPM server.'); end
n       = regexp(s,'spm8_updates_r(\d.*?)\.zip','tokens','once');
if isempty(n)
    fprintf('         There are no updates available yet.\n');
    return;
else
    n   = str2double(n{1});
end

try
    [v,r] = spm('Ver','',1); r = str2double(r);
catch
    error('SPM cannot be found in MATLAB path.');
end
if ~strcmp(v,'SPM8'), error('Your SPM version is %s and not SPM8',v); end
rs = [3042 3164 3408 3684 4010 4290 4667 5236 6313];
if isnan(r), r = rs(1); end 
if floor(r) == 8
    try
        r = rs(round((r-floor(r))*10)+1);
    catch
        r = rs(end);
    end
end

if n > r
    fprintf('         A new version of SPM is available on:\n');
    fprintf('     %s\n',url);
    fprintf('        (Your version: %d - New version: %d)\n',r,n);

    if update
        d = spm('Dir'); 
        delete(get(0,'Children')); spm('clean'); evalc('spm_rmpath');
        try
            s = unzip([url sprintf('spm8_updates_r%d.zip',n)], d);
            fprintf('             %d files have been updated.\n',numel(s));
        catch
            fprintf('          Update failed: check file permissions.\n');
        end
        addpath(d);
    end
else
    fprintf('         Your version of SPM is up to date.\n');
end
