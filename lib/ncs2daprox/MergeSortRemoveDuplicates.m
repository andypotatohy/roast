function [fx,fy,index]=MergeSortRemoveDuplicates(fxSet1,fySet1,Set1index,fxSet2,fySet2,Set2index)

% For two sets of feature points: Merge, sort and remove duplicates 
% Sorting is based on indices of mergerd set.


% All input/output paramters (i.e. fxSet1,fySet1,Set1index,fxSet2,fySet2,Set2index,fx,fy,index) should be
% row vectors or column vector.
% For example fxSet1 in column vector format
% x1
% x2
% ...
% xn

% For example fxSet1 in row vector format
% x1,x2,...,xn


[r1 c1]=size(fxSet1);
[r2 c2]=size(fySet1);
[r3 c3]=size(Set1index);
[r4 c4]=size(fxSet2);
[r5 c5]=size(fySet2);
[r6 c6]=size(Set2index);


if (r1==1)
    if( r1~=r2 ||  r1~=r3 ||  r1~=r4 || r1~=r5 || r1~=r6 )
    disp('Message from MergeSortRemoveDuplicates.m');
    disp('Invalid Data, all data must be either row vector or column vector');
    return
    end
elseif (c1==1)
    if( c1~=c2 ||  c1~=c3 ||  c1~=c4 || c1~=c5 || c1~=c6  )
    disp('Message from MergeSortRemoveDuplicates.m');
    disp('Invalid Data, all data must be either row vector or column vector');
    return
    end
else
    disp('Message from MergeSortRemoveDuplicates.m');
    disp('Invalid Data, some data is not vector');
    return
end



if (r1==1) % row vector
    A=[Set1index' fxSet1' fySet1']; 
    B=[Set2index' fxSet2' fySet2'];
    
else % column vector
    A=[Set1index fxSet1 fySet1]; 
    B=[Set2index fxSet2 fySet2];
end

% A and B has following format 
%[index1 x1 y1;
% index2 x2 y2;
%   ,...,     ; 
% indexn xn yn]

% c = union(A,B,'rows') when A and B are matrices with the same number of
% columns returns the combined rows from A and B with no repetitions.

C = union(A,B,'rows');
C = sortrows(C,1); % sorts the matrix based on the columns specified (here 1)

index=C(:,1); % move back
fx=C(:,2);
fy=C(:,3);

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------

