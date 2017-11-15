close all, clc, clear all 

%%%%%%%%%%%%%%%%%%%%%%%%% SINE APPROXIMATION %%%%%%%%%%%%%%%%%%%%%%
disp('Sine Approximation using Natural Cubic Spline');
x = -pi:.05:pi; x=x';
y = sin(x);   
size(x)
ei= length(x);
initialbreaks=[1;ei]; % user provided indices of initial break points (or corer points/control points)
%initialbreaks=[1;floor(ei/2);ei]; % OR first, middle, last
MaxAllowedDist=.05; % Max. allowed Distance b/w original data and fitted parametric spline
MaxAllowedSqDist=MaxAllowedDist^2 % Max. allowed Square Distance... 


%%%  Some possible combinations to call the function ncs2dapprox.
%%  Use only one , comments others


%[bx,by,finalbreaks]=ncs2dapprox(x,y)                             % 1     

%[bx,by,finalbreaks]=ncs2dapprox(x,y,MaxAllowedSqDist)            %or 2  
 
%[bx,by, finalbreaks]=ncs2dapprox(x,y,MaxAllowedSqDist)           %or 3  (recommended call)   

%[bx,by,finalbreaks]=ncs2dapprox(x,y,MaxAllowedSqDist,initialbreaks)      % or 4    

%[bx,by,finalbreaks]=ncs2dapprox(x,y,MaxAllowedSqDist,initialbreaks)       % or 5

[bx,by,finalbreaks,MasSqDist]=ncs2dapprox(x,y,MaxAllowedSqDist)  % or 6 (recommended call)   

 %[bx,by,finalbreaks,MasSqDist]=ncs2dapprox(x,y,MaxAllowedSqDist,initialbreaks)  % or 7


%%% Now we can only save break points/control points i.e. bx,by and their indices i.e. finalbreaks, optionally  MasSqDist and any time we can
%%% generate the approximation of original data by spline interpolation
%%% using finalbreaks,bx,by

%%% Now we would use the control points returned by fuction to generate
%%% approximationg Spline using spline interpolation

t = finalbreaks';    
pp1= spline(t,[bx,by]');
range = linspace(1,finalbreaks(end),finalbreaks(end)); 
xiyi = ppval(pp1,range);
xi=xiyi(1,:)';
yi=xiyi(2,:)';
    
figure ,hold on
plot(bx,by,'ro','LineWidth',2)
plot(x,y,'b','linewidth',2)
plot(xi,yi,'g','linewidth',2)
title('\bf{<Blue: Original Data>,<Green: Spline Approximation>,<Red Circle: Break Point>}');

%%%%%%%%%%%%%%%%%%%%% CIRCLE APPROXIMATION %%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------------------------');
disp('Circle Approximation using Natural Cubic Spline');
 
rad=1;
theta = linspace(0, 2*pi, 300);
x = rad*cos(theta);
y = rad*sin(theta);

ei= length(x);
initialbreaks=[1;ei];
MaxAllowedDist=.05; % Max. allowed Distance b/w original data and fitted parametric spline
MaxAllowedSqDist=MaxAllowedDist^2 % Max. allowed Square Distance... 

[bx,by,finalbreaks,MasSqDist]=ncs2dapprox(x,y,MaxAllowedSqDist,initialbreaks)      

t = finalbreaks';    
pp1= spline(t,[bx,by]');
range = linspace(1,finalbreaks(end),finalbreaks(end)); 
xiyi = ppval(pp1,range);
xi=xiyi(1,:)';
yi=xiyi(2,:)';

figure ,hold on
plot(bx,by,'ro','LineWidth',2)
plot(x,y,'b','linewidth',2)
plot(xi,yi,'g','linewidth',2)
axis equal
title('\bf{<Blue: Original Data>,<Green: Spline Approximation>,<Red Circle: Break Point>}');

%%%%%%%%%%%%%%%%%%%%%%%%%  FIVE APPROXIMATION  %%%%%%%%%%%%%%%%%%%%%%%
disp('----------------------------------------------');
disp('Five Approximation using Natural Cubic Spline');

M = dlmread('five.txt'); % read from ASCII file
x=M(:,1);
y=M(:,2);
clear M
ei= length(x);
initialbreaks=[1;ei];
MaxAllowedDist=4; % Max. allowed Distance b/w original data and fitted parametric  spline
MaxAllowedSqDist=MaxAllowedDist^2 % Max. allowed Square Distance... 

[bx,by,finalbreaks,MasSqDist]=ncs2dapprox(x,y,MaxAllowedSqDist)     


t = finalbreaks';        % t        is 1 x n column vector
pp1= spline(t,[bx,by]'); % [bx,by]' is 2 x n column vector
[breaks,coefs,l,order,d] = unmkpp(pp1);
range = linspace(1,finalbreaks(end),finalbreaks(end)); 
xiyi = ppval(pp1,range);
xi=xiyi(1,:)';
yi=xiyi(2,:)';

figure ,hold on
plot(bx,by,'ro','LineWidth',2)
plot(x,y,'b','linewidth',2)
plot(xi,yi,'g','linewidth',2)
title('\bf{<Blue: Original Data>,<Green: Spline Approximation>,<Red Circle: Break Point>}');


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------