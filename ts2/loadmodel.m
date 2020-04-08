function loadmodel(modelname);

%function loadmodel(modelname)
%
%load a radial basis model from disk
%
% M. Small 
% Created: 7/6/99
% Updated: 7/6/99

rb_get_globals

if nargin<1,
   modelname='model';
end;

eval(['load ',modelname]);
disp(['loading : ',modelname]);
