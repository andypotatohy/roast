Approximation of 2 Dimension Data by Natural Cubic Spline
-------------------------------------------------------
(Impatient user Run the program Testncs2dapprox.m)

Suppose we have set of continuous points (xi,yi), 1<=i<=n (e.g. boundary or some signal) and we want to approximate them using Natural Cubic Spline. This parametric representation of data by spline has many advantages like compaction, smoothness etc.

A general concept of fitting Algorithm is following:

1.Fit the spline to Data using initial break points.
2.Find the Max. square distance b/w spline and data.
3.while(Max. Square Distance b/w Spline & Data to be fitted > Max Allowed Square Distance)    
4.  Add point of max. distance to set of break points.
5.  Fit the spline using new set of break points.
6.  Find the Max. square distance b/w spline and data.  
7.end


Two files are most important and brief description about them is following.

ncs2dapprox.m
-------------
Core function that do spline approximation of data as explained in above algorithm using uniform parameterization. 
This function can be called with variable number of input/output arguments.Call with at least two input arguments and three output arguments.

Syntax of Usage:
[arg1out,arg2out,arg3out,arg4out]=ncs2dapprox(arg1in,arg2in,arg3in,arg4in)


arg1in: Input x-Data e.g. [x1, x2, x3,...,xn]
arg2in: Input y-Data e.g. [y1, y2, y3,...,yn]
arg3in: Maximum allowed Square Distance between Data and parametric values (Optional argument)
arg4in: Indices of Data where Spline MUST interpolate (Optional argument)

arg1out: x-values of output break points
arg2out: y-values of output break points
arg3out: Indices of output break points
arg4out: max squared distance b/w input and output values

See Testncs2dapprox.m how we call this function using various options.


Testncs2dapprox.m
-----------------
A Test program that shows how to use ncs2dapprox.m , Some of the possible calls are
commented, you can uncomment and use them.

% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % --------------------------------


