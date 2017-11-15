function [squaredmax,rowIndex]=MaxSqDistAndRowIndexbw2Mat(mat1,mat2)
% mat1=[x11,x12,...x1n;
%       x11,x12,...x1n; 
%       ...
%       xn1,xn2,...xnn; 

% mat2=[y11,y12,...y1n;
%       y11,y12,...y1n; 
%       ...
%       yn1,yn2,...ynn; 
% OR in brief mat1 and mat2 format is like following
%                               [P1;
%                                P2;
%                                P3;
%                                P4;
%                                ...
%                                PN];

%% we have to find square distance between each row of mat1 and mat2 and then among
%% all distances we have to return maximum one and its row index in mat1 (or
%% mat2, same for both)

% returns -1 for squaredmax and rowIndex for exception cases
% the algorithms is based on euclidean distance 

[r1 c1]=size(mat1);
[r2 c2]=size(mat2);
squaredmax=-1;
rowIndex=-1;

%%% Exception Check
if ( r1~=r2 || c1~=c2  )
    disp('Message from MaxSqDistAndRowIndexbw2Mat.m');
    disp('Two matrices must be of equal size');
    return
end

%%% empty matrices check
if ( r1==0  )
    disp('Message from MaxSqDistAndRowIndexbw2Mat.m');
    disp('Empty matrices');
    return
end

%%%%% Actual Algorithm
for k=1:r1
    SqDist=TwoNormSqDist(mat1(k,:),mat2(k,:)); % square distance b/w kth row
    if(SqDist > squaredmax )
        squaredmax=SqDist;
        rowIndex=k;
    end
end


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------




