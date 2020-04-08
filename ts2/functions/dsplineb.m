function y=dsplineb(x,p);

% Brian Fleming, 11/6/00.

global SPLINE_ORDER
shift=(SPLINE_ORDER+1)/2;

y=dbspline(x,p-shift);
