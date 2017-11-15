% return 1 if input argument is vecotor else return 0

function ans=isvec(x)
ans=0;
d=size(x);

if(length(d)>2) % not vector
    return 
end

[r c]=size(x);

if (r>1 & c>1) % not vector
    return
end

ans=1; % vector


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------
    

