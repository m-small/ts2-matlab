function rb_viewbase(base,v,func);

%function rb_viewbase(base,v,func);
%
% show parameters associated with a set of basis functions.
%
%
%
% M. Small 
% Created: 3/4/98
% Updated: 20/8/98

disp(' ');
disp(base)

disp(' ');
[nc,d]=size(base.centres);
[nr,d]=size(base.radii);
numbers=num2str([base.centres;base.radii],4);
disp('base.centres =');
disp([ones(nc,1)*'    ',numbers(1:nc,:)]);
disp(' ');
[nr,d]=size(base.radii);
disp('base.radii =');
disp([ones(nr,1)*'    ',numbers((nc+1):(nc+nr),:)]);
disp(' ');
disp('v(base.strategy),:)'' =');
v=v(base.strategy,:)';
[nv,d]=size(v);
for i=1:nv,
   disp(sprintf('      %3d ',v(i,:)));
end;
disp(' ');
disp('func{base.func(:)} =');
str=[];
for i=1:d,
   thisone=[blanks(9),func{base.func(i)}];
   thisone=thisone((end-8):end);
   str=[str,thisone,' '];
end;
disp(str);
