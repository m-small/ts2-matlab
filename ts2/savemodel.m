function savemodel(modelname);

%function savemodel(modelname)
%
%save a radial basis model to disk
%
% M. Small 
% Created: 7/6/99
% Updated: 7/6/99

rb_get_globals

if nargin<1,
   modelname='model';
end;

eval(['save ',modelname,' rb_method rb_x rb_y rb_base rb_lambda rb_delta rb_basis rb_penalty rb_descr_length rb_error rb_embed rb_functions rb_timescale']);
disp(['saving : ',modelname]);
