function base=rb_catbase(base1,base2,base3);

%
% M. Small 
% Created: 5/4/98
% Updated: 5/4/98

if nargin<3,
   base3=rb_nullbase;
end;

base.centres=[base1.centres base2.centres base3.centres];
base.strategy=[base1.strategy; base2.strategy; base3.strategy;];
base.func=[base1.func; base2.func; base3.func;];
base.radii=[base1.radii base2.radii base3.radii];

