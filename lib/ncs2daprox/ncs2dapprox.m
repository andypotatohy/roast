
function [bx,by,finalbreaks,MaxSqDist]=ncs2dapprox(x,y,varargin)

s1='Function Calling Syntax. Call with atleast two input arguments and three ouput arguments';
s2='[arg1out,arg2out,arg3out,arg4out]=ncs2dapprox(arg1in,arg2in,arg3in,arg4in)';
s3='arg1in: Input Data e.g. [x1, x2, x3,...,xn]';
s4='arg2in: Input Data e.g. [y1, y2, y3,...,yn]';
s5='arg3in: Maximum allowed Square Distance between Data and paramteric values (Optional argument)';
s6='arg4in: Indices of Data where Spline MUST interpolate (Optional argument)';
s7='arg1out: x-values of output break points';
s8='arg2out: y-values of output break points';
s9='arg3out: Indices of ouptput break points';
s10='arg4out: max squared distance b/w input and output values';
if(nargin<2 || nargin>4 ) 
    disp(s1);   disp(s2);   disp(s3);    disp(s4);    disp(s5);
    disp(s6);   disp(s7);   disp(s8);    disp(s9);    disp(s10);
    return
end   
if(nargout<3 ) 
    disp(s1);   disp(s2);   disp(s3);    disp(s4);    disp(s5);
    disp(s6);   disp(s7);   disp(s8);    disp(s9);    disp(s10);
end

if(~isvec(x) || ~isvec(y) )
    error('first two input arguments must be row OR column vector');    
end
[x]=getcolvector(x);     % make column vector, if it is row vector
[y]=getcolvector(y);
if(length(x)~=length(y))
    error('first two arguments must have equal number of values');    
end

if (size(x,1)<2)
    error('Atleast two values are required in data (first two arguments)');    
end

%%% Default Values %%%
MaxAllowedSqDist=1;
initialbreaks=[1; size(x,1)]; % first & last row indices 
defaultValues = {MaxAllowedSqDist initialbreaks};
%%% Default Values %%%

%%% Assign Values %%%
nonemptyIdx = ~cellfun('isempty',varargin);
defaultValues(nonemptyIdx) = varargin(nonemptyIdx);
[MaxAllowedSqDist initialbreaks] = deal(defaultValues{:});
%%% Assign Values %%%

finalbreaks=initialbreaks; % initially same

finalbreaks=getcolvector(finalbreaks); % make column vector, if it is row vector

bx=x(finalbreaks);
by=y(finalbreaks);

% appending first,last data values and indices in
% bx,by,finalbreaks
[bx,by,finalbreaks]=MergeSortRemoveDuplicates(bx,by,finalbreaks,x(1),y(1),1); 
[bx,by,finalbreaks]=MergeSortRemoveDuplicates(bx,by,finalbreaks,x(end),y(end),size(x,1)); 

t = finalbreaks';
pp1= spline(t,[bx,by]');

range = linspace(1,length(x),length(x)); % uniform parameterizaton 
xiyi = ppval(pp1,range);

xi=xiyi(1,:)';
yi=xiyi(2,:)';

[MaxSqDist, MaxSqDistIndex]=MaxSqDistAndRowIndexbw2Mat([x,y],[xi,yi]); % estimate error

while(MaxSqDist > MaxAllowedSqDist)    
    [bx,by,finalbreaks]=MergeSortRemoveDuplicates(bx,by,finalbreaks,x(MaxSqDistIndex),y(MaxSqDistIndex),MaxSqDistIndex);

    t = finalbreaks';    
    pp1= spline(t,[bx,by]');
    
    range = linspace(1,length(x),length(x));
    xiyi = ppval(pp1,range);
    
    xi=xiyi(1,:)';
    yi=xiyi(2,:)';
    
    [MaxSqDist, MaxSqDistIndex]=MaxSqDistAndRowIndexbw2Mat([x,y],[xi,yi]);
    
end





% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------

