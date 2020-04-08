function y=splineb(x,p);

% Brian Fleming, 11/6/00.
  
global SPLINE_ORDER
shift=SPLINE_ORDER/2;
p(2,:)=p(2,:)-shift*p(1,:);

y=bspline(x,p);

