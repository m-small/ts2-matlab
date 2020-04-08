function rb_viewmodel

%function rb_viewmodel
%
%print out the parameters of the current radial basis model (stored as global variables).
%
%M. Small 
% Created: 6/6/99
% Updated: 6/6/99

rb_get_globals

disp('size(rb_x)=');
disp(size(rb_x));
disp('size(rb_y)=');
disp(size(rb_y));
disp('rb_penalty=');
disp(rb_penalty);
disp('rb_descr_length=');
disp(rb_descr_length);
disp('rms(rb_error)=');
disp(rms(rb_error'));
lambdas(rb_basis,:)=rb_lambda;
deltas(rb_basis)=rb_delta;
disp('[rb_lambda rb_delta (scales)]=');
if ~isempty(lambdas),
   if length(lambdas(1,:))>1,
      disp([lambdas(:,1) deltas' lambdas(:,2)]);
   else,
      disp([lambdas deltas']);
   end;
end;
disp('rb_base=');
rb_viewbase(rb_base,rb_embed,rb_functions);

