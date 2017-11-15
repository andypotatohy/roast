% if vin is row vector change it to column vector then return it.
% if vin is column vector then return it as it is.
% if vin is not a vector diplay error message and  returns control to the
% keyboard
function [vout]=getcolvector(vin)

if(~isvec(vin)  )
    error('getrowvector.m => input must be a vector');    
end

vout=vin;
[r c]=size(vin);
if (r==1)      % make row column if it is row vector
    vout=vin';
end


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------
