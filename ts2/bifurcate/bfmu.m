function mu=bfmu(mu);

%rescale bifurcation parameter values using global model.
%
% mu=[] or omitted, then mu is a 2-by-1 vector of initial value and increment.
% otherwise mu (out) is mu (in) rescaled according to rb_x
% initial value and increment are deduced from the last row of
% rb_x, if they are not linearly scaled the initial value and
% increment are given that for the number of points in rb_x would
% range from minimum to maximum value in equal steps.
%
% if mu=-1 then mu (output) is the values stored as the last row of 
% rb_x. 
% 
%
% M. Small 
% Created: 18/10/99
% Updated: 18/10/99

if nargin<1,
   mu=[];
end;

rb_get_globals;

[dx,nx]=size(rb_x);

minmu=min(rb_x(dx,:));
maxmu=max(rb_x(dx,:));
nmu=nx;

if isempty(mu)
   mu(1)=minmu;
   mu(2)=(maxmu-minmu)/(nmu-1);
else,
  if max(size(mu))==1 & mu(1)==-1,
     mu=rb_x(dx,:);
else,  
   mu=minmu+mu*((maxmu-minmu)/(nmu-1));
  end;
end;

