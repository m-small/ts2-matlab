%script file rb_clear_globals
%
% rb_get_globals gets the global variables associated with a rb model.
%
% M. Small 
% Created: 6/6/98
% Updated: 6/6/98

%global pl_method 					% method used in model

%the data the model is fitted to
clear global rb_method
clear global rb_x 							% independent variables
clear global rb_y 							% dependent variables.
clear global rb_base
clear global rb_lamba 	          		% The weights
clear global rb_delta						% The (absolute) precison
clear global rb_basis
clear global rb_penalty
clear global rb_descr_length 				% The actual description length of this model
clear global rb_error 						% y - rb_phi(rb_base) * lambda
clear global rb_embed						% the embedding strategies
clear global rb_functions
clear global rb_timescale 					% timescale used in pl_spacetime
