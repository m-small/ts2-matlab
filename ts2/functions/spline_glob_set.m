function spline_glob_set(order)
% spline_glob_set(order)
%
% Sets up global variables for spline functions (default order=3)
%
% We set up the global variables 'SPLINE_ORDER' and 'SPLINEWAVE_COEFFS'. These
% are set for the following functions:
%
%			bspline.f (mex)
%		       dbspline.m
%			bsplinewave.f (mex)
%		       dbsplinewave.f (mex)
%
% Brian Fleming, 14/6/00.

if nargin<1,
  order = 3;
end;


global SPLINE_ORDER
SPLINE_ORDER = order;

global SPLINEWAVE_COEFFS
SPLINEWAVE_COEFFS = bwave_coeffs(SPLINE_ORDER);
