% if vin is column vector change it to row vector then return it.
% if vin is row vector then return it as it is.
% if vin is not a vector diplay error message and  returns control to the
% keyboard
function [vout]=getrowvector(vin)

if(~isvector(vin)  )
    error('getrowvector.m => input must be a vector');    
end

vout=vin;
[r c]=size(vin);
if (c==1)      % make row vector if it is column vector
    vout=vin';
end

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------

