function y=waveb(x,p);

% Brian Fleming, 11/6/00.

global SPLINE_ORDER
shift=SPLINE_ORDER-0.5;
p(2,:)=p(2,:)-shift*p(1,:);

y=bsplinewave(x,p);
