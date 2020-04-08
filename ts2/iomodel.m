function iomodel(x_fit,y_fit,td,func,penalty,niterat,nbuilds,nbasis,nexpand,wobble,meth);

%function iomodel(x,y,testdatum,func,penalty,niterat,nbuilds,nbasis,nexpand,wobble,method);
%
% build a model f such that f(x)-y is iid. model is built
%
% func is the basis functions to use. current choices include;
%  'cubic','quintic','gaussian','tophat','wavelet' and 'morlet'
% penalty is a string specifying the penalty equation to eval e.g.
%   Akaike = 'dim_x*log(mss)+2*k'                                          --- Akaike criterion
%   Schwarz = 'dim_x*log(mss)+k*log(dim_x)',                               --- Schwarz criterion
%   Rissanen = 'dl1(mss,a,deltas,dim_x)';                                  --- pseudo-linear DL
%   Rissanen2= 'dl2(mss,a,deltas,dim_x,base,basic,offset,v)';              --- approx DL
% method is a string composed of either 
%    'r' (radial basis functions),
%    'w' (optimal wavelet filter), OR
%    'n' (1-layer neural network)  AND
% none, some or all of 'c' (constant term) and 'l' (linear terms)
% if wobble=0 then no nonlinear optimisation of selected bases is
% performed, if wobble>0 then it is and wobble is a parameter of
% the penalty function that prevents basis functions from drifting
% far from the data. if models are crap, lower wobble (wobble=0.1
% is a good first guess).
%
% Default
%  func={'cubic','quintic','sigmoid','gaussian','tophat','wavelet','morlet'}
%  penalty=Rissanen
%  niterat=2
%  nbuilds=3
%  nbasis=100
%  nexpand=5
%  wobble=0
%  method='clr'
%
% M. Small 
% Created: 6/6/99
% Updated: 15/2/02


%possible penalties
Akaike = 'dim_x*log(mss)+2*k';                %AIC
Schwarz = 'dim_x*log(mss)+k*log(dim_x)';      %SIC/BIC
%Rissanen = 'dl1(mss,unscaled_a,unscaled_deltas,dim_x)';  %pseudo-linear DL
Rissanen = 'dl1(mss,a,deltas,dim_x)';  %pseudo-linear DL
%i.e. penalise for description of weights only
%Rissanen2= 'dl2(mss,unscaled_a,unscaled_deltas,dim_x,base,basic,offset,v)'; %approx DL
Rissanen2= 'dl2(mss,a,deltas,dim_x,base,basic,offset,v)'; %approx DL
%assume penalty for all significant parameters are equal and equal to the penalty for
%the weights (this assumption is not necessarily well founded --- just convenient)



%check for inputs
na=nargin;
if na<11,
  meth=[];
end;
if na<10,
   wobble=[];
end;
if na<9,
   nexpand=[];
end;
if na<8,
   nbasis=[];
end;
if na<7,
   nbuilds=[];
end;
if na<6,
   niterat=[];
end;
if na<5,
   penalty=[];
end;
if na<4,
   func=[];
end;
if na<3,
   td=[];
end;
if na<2,
   y_fit=[];
end;
if na<1,
   x_fit=[];
end;
clear na

%default values
if isempty(meth),
  meth='clr';
  disp(['method=''',meth,''';']);
end;
if isempty(wobble),
   wobble=1;
   disp(['wobble=',int2str(wobble),';']);
end;
if isempty(nexpand),
   nexpand=5;
   disp(['nexpand=',int2str(nexpand),';']);
end;
if isempty(nbasis),
   nbasis=100;
   disp(['nbasis=',int2str(nbasis),';']);
end;
if isempty(nbuilds),
   nbuilds=3;
   disp(['nbuilds=',int2str(nbuilds),';']);
end;
if isempty(niterat),
   niterat=2;
   disp(['niterat=',int2str(niterat),';']);
end;
if isempty(penalty),
   penalty=Rissanen2;
   disp(['penalty=''',penalty,''';']);
end;
if isempty(func),
   func={'cubic','quintic','sigmoid','gaussian','tophat','wavelet','morlet'};
   disp('func=');
   disp(func);
end

%embed data, sort out y_test and y_fit
[dx,nx]=size(x_fit);
y_fit=y_fit(:)';
[dy,ny]=size(y_fit);
if ny~=nx,
   error('size of x_fit and y_fit don''t match');
end;
if dx<15,
	v=genve([-1 0:(dx-1)]);
else,
   v=[-1 0:(dx-1)];
end;
if isempty(td),
   td=floor(nx/2);
end;
if td>nx,
   error('testdatum too big');
end;
if ~isnan(td),
	y_test=y_fit((td+1):end);
	x_test=x_fit(:,(td+1):end);
	y_fit=y_fit(1:td);
	x_fit=x_fit(:,1:td);
end;


%initialise model parameters
best_mdl=inf;
best_base=[];
best_lambda=[];
best_delta=[];
best_err=[];

%check for current model
pbase=[];mdl=inf;
rb_get_globals
samemeth=1;
for i=1:length(rb_method),
  samemeth=samemeth & any(meth==rb_method(i));
end;
if samemeth & length(meth)==length(rb_method),
  if all(size(rb_functions)==size(func)),
    if all(size(rb_embed)==size(v)),
      if (rb_embed(~isnan(rb_embed(:)))==v(~isnan(v(:)))),
	disp('Using basis functions from current global model as candidates.');
	pbase=rb_base;
	mdl=rb_descr_length;
	improved=0;
	%present model is the best model so far ---  for now.
	best_mdl=mdl;
	best_base=rb_base;
	best_lambda=rb_lambda;
	best_delta=rb_delta;
	best_basis=rb_basis;
	best_err=rb_error;
if ~isempty(best_lambda),
	disp(['* Current best model : ',int2str(length(best_lambda(:,1))),' basis functions, RMS =',num2str(rms(best_err')),', MDL =',num2str(best_mdl),'.']);
end;
end;
    end;
  end;
end;
clear samemeth

%build the model
for i=1:niterat
   base=pbase;
   for j=1:nbuilds,
      if trace,
         disp(['Starting calculation for node (',int2str(i),',',int2str(j),')']);
      end;
      [base,lambda,delta,basis,err,mdl]=rb_topdown(x_fit,y_fit,v,func,base,penalty,nbasis,nexpand,wobble,meth);
      disp(['Calculation complete for node (',int2str(i),',',int2str(j),')']);
      %is latest model the best?
      if mdl<best_mdl
         improved=1;
         if isempty(lambda),
             disp(['* Current best model : 0 basis functions (null model), RMS =',num2str(rms(err')),', MDL =',num2str(mdl),'.']);
        else,
             disp(['* Current best model : ',int2str(length(lambda(:,1))),' basis functions, RMS =',num2str(rms(err')),', MDL =',num2str(mdl),'.']);
         end;
         best_mdl=mdl;
         best_base=base;
         best_lambda=lambda;
         best_delta=delta;
         best_basis=basis;
         best_err=err;
         if trace,
            ylabel(['(',int2str(i),',',int2str(j),') of (',int2str(niterat),',',int2str(nbuilds),')']);
         end;
         %and store best model in globals
         rb_set_globals(x_fit,y_fit,best_base,best_lambda,best_delta,best_basis,penalty,best_mdl,best_err,v,func,meth);
      elseif trace,
         close;
      end;
   end;
end;

if ~improved,
   disp('Didn''t improve initial model --- bummer.');
end;



if isnan(td),
   %do nothing.
elseif ~isempty(y_test),
   %one step predictions
   figure;
   %dishonest predictions
yp=rb_eval(x_fit)';
ep=yp-y_fit;
yt=y_fit;   
subplot(211);
plotcols(ep,'b:',yt,'y-',yp,'r--');
zoom on;
title(['dishonest predictions RMS = ',num2str(rms(ep)),' (STD = ',num2str(std(yt)),')']);
   %honest predictions

yp=rb_eval(x_test)';
ep=yp-y_test;
yt=y_test;   
subplot(212);
plotcols(ep,'b:',yt,'y-',yp,'r--');
zoom on;
title(['honest predictions RMS = ',num2str(rms(ep)),' (STD = ',num2str(std(yt)),')']);
xlabel('true (yellow), predicted (red), and error (blue)');
   %honest predictions
 end;
