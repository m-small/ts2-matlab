function [a,rang]=rb_param(func);

%function [a,rang]=rb_param(function);
%
% returns the number and range of the parameters associated with 
% a particular basis function.
% 
% current, valid basis functions include;
%  gaussian
%  tophat
%  sigmoid
%  cubic 
%  quintic
%
% if one bound is nan then the values are gaussian with standard
% deviation of the other bound (time the spread of the data), if
% both are nan then they have standard deviation of the spread of
% the data. Otherwise candidate radii are uniformly selected
% between the two bounds.
%
% M. Small 
% Created: 31/3/98
% Updated: 7/7/00

if ~isstr(func),
   disp('WARNING: non string data type in rb_param');
end;

switch lower(func),
case 'gaussian', a=1;rang=[2;nan;];
case 'wavelet', a=1;rang=[2;nan;];
case 'morlet', a=2; rang=[[2 4];[nan 6];];
case 'tophat', a=2;rang=[[2 1];[nan 6];];
case 'sigmoid', a=1;rang=[1;nan];
%case 'sigmoid', a=2;rang=[[0 -pi];[3 0]];
case 'cubic', a=0;rang=[];
 case 'quintic', a=0;rang=[];
 case 'bspline', a=2;rang=[[eps 0];[2 0];];       %big splines for trends
 case 'bsplinewave', a=2;rang=[[eps 0];[0.1 0];]; %lil' wavelets for wiggles
 case 'splineb', a=2;rang=[[eps 0];[2 0];]; % ditto
 case 'waveb', a=2;rang=[[eps 0];[0.1 0];]; %  " "
    case 'doublet', global TP_1D_DIV; a=1;rang=[eps; 10/TP_1D_DIV];%a=2;rang=[[2 nan;];[nan nan;];];
    case 'singlet', global TP_1D_DIV; a=1;rang=[eps; 10/TP_1D_DIV];%a=2;rang=[[2 nan;];[nan nan;];];
otherwise, disp(['WARNING: invalid function type - ',func,' - in rb_param']);
end;

