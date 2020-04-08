function rb_view2d(v);

%function rb_view2d(v)
%
%2-D graphic representation of the basis functions of a model
%
%
%
% M. Small 
% Created: 20/11/01
% Updated: 20/11/01

if nargin<1,
  v=[0 1];
end;
v=v+1;

rb_get_globals;
plotmat(rb_x(v,:),'k.');
hold on;
nb=length(rb_base.func);

if any(rb_method=='r'),
    plotmat(rb_base.centres(v,:),'r*');
    xt=cos(0:0.1:6.3);
    yt=sin(0:0.1:6.3);
    for i=1:nb,
        plot(rb_base.centres(v(1),i)+rb_base.radii(i)*xt,...
            rb_base.centres(v(2),i)+rb_base.radii(i)*yt,'r');
    end;
elseif any(rb_method=='n'),
    for i=1:nb,
        plot([rb_base.radii(1,i)/rb_base.centres(v(1),i) 0],...
            [0 rb_base.radii(1,i)/rb_base.centres(v(2),i)],'r');
     end;
 end;
 
