% This is a mex file for MATLAB. It calculates the values of the cardinal 
% B-spline scaling function for a given (vector) set of parameters (for a 
% vector of input values), used to specify the dilation, translation and 
% order parameter of the function. The function is normalised to have 
% constant L^2 norm for different dilation and translation paramters.
%  See "Fundamentals of Wavelets", pp.106 by Goswami abd Chan.
%
%  INPUT:
%   'x' vector of input arguments to function
%   'pms' vector of function parameters	
%
%  OUTPUT:
%   'y' function values
%
%  Notes:
%   'x' is a 1-by-n or n-by-1 vector
%   'y' is 1-by-n
%    y(i) is the evaluation of the function at x(i) with parameters values
%    pms(:,i) (or pms(:) if pms is 2-by-1).
%
% The format of the vector 'pms' is as follows:
%    pms = [a; b];
% with
%   'a' dilation parameters (size = [1 n])
%   'b' translation parameters (size = [1 n])
%
% The order of the spline, parameter 'ord', is taken in from a global 
% variable 'SPLINE_ORDER' in the global MATLAB workspace.
%   'ord' order of spline (polynomial of degree ord-1) 
%
% NOTE: If a vector 'x' is given with only one paramter 'ord'. We 
% calculate the output of the B-spline for values 'x' and order 'ord'.
%
% Brian Fleming, 11/6/00.
