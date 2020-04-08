function rb_set_globals(x,y,base,lambda,delta,basis,penalty,mdl,err,v,func,meth);

%function rb_set_globals(x,y,base,lambda,delta,basis,penalty,mdl,err,v,func,meth);
%
% rb_globals gets the global variables associated with a rb model.
%
% M. Small 
% Created: 6/6/99
% Updated: 6/6/99

rb_get_globals
rb_x=x;
rb_y=y;
rb_base=base;
rb_lambda=lambda;
rb_delta=delta;
rb_basis=basis;
rb_penalty=penalty;
rb_descr_length=mdl;
rb_error=err;
rb_embed=v;
rb_functions=func;
rb_method=meth;
