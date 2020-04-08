function y=dwaveb(x,p);
  
global SPLINE_ORDER
shift=SPLINE_ORDER-0.5;
  
y=dbsplinewave(x,p-shift);
