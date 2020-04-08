function base=rb_assbase(base,i);

%
% M. Small 
% Created: 5/4/98
% Updated: 5/4/98

base.centres=base.centres(:,i);
base.strategy=base.strategy(i);
base.func=base.func(i);
base.radii=base.radii(:,i);