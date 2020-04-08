%script file rb_get_globals
%
% rb_get_globals gets the global variables associated with a rb model.
%
% M. Small 
% Created: 30/3/98
% Updated: 6/6/99

global rb_method 					% method used in model
if isempty(rb_method),
  rb_method='clr';
end; %backwards compatibility


%the data the model is fitted to
global rb_x 							% independent variables
global rb_y 							% dependent variables.

%the basis functions
global rb_base
%global rb_base.centres 			% Coordinates of the centre
%global rb_base.radii 		 		% Radii of centres
%global rb_base.strategy 			% Embedding strategy at centre
%global rb_base.func 				% The basis function used

%model weights and precisions
global rb_lambda 	          		% The weights
global rb_delta						% The (absolute) precison
global rb_basis                  % The relevant columns


%the models description length
global rb_penalty
global rb_descr_length 				% The actual description length of this model

%the models error
global rb_error 						% y - rb_phi(rb_base) * lambda

%the embedding(s)
global rb_embed						% the embedding strategies

%the basis function(s)
global rb_functions

%timescale?
global rb_timescale 					% timescale used in pl_spacetime

%global pl_ellipse               % global variable for ellipsoid basis functions.
