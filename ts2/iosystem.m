function iosystem(x_fit,y_fit,x_test,y_test,v,func,penalty,niterat,nbuilds,nbasis,nexpand,wobble,meth);

%function iosystem(x_fit,y_fit,x_test,y_test,v,func,penalty,niterat,nbuilds,nbasis,nexpand,wobble,method);
%
% func is the basis functions to use. current choices include;
%  'cubic','quintic','gaussian','tophat','wavelet' and 'morlet'
% penalty is a string specifying the penalty equation to eval e.g.
%   Akaike = 'dim_x*log(mss)+2*k'                                          --- Akaike criterion
%   Schwarz = 'dim_x*log(mss)+k*log(dim_x)',                               --- Schwarz criterion
%   Rissanen = 'dl1(mss,a,deltas,dim_x)';                                  --- pseudo-linear DL
%   Rissanen2= 'dl2(mss,a,deltas,dim_x,base,basic,offset,v)';              --- approx DL
%
% Default
%  v=[-1 0 1 2]
%  func={'cubic','quintic','sigmoid','gaussian','tophat','wavelet','morlet'}
%  penalty=Rissanen
%  niterat=2
%  nbuilds=3
%  nbasis=100
%  nexpand=5
%  wobble=0
%  method='clr'
%
% designed to work for input/output type systems, (i.e. not time
% series, not embedded). Not tested.
%
% M. Small 
% Created: 6/6/99
% Updated: 5/7/00


%possible penalties
Akaike = 'dim_x*log(mss)+2*k';                %AIC
Schwarz = 'dim_x*log(mss)+k*log(dim_x)';      %SIC/BIC
%Rissanen = 'dl1(mss,unscaled_a,unscaled_deltas,dim_x)';  %pseudo-linear DL
Rissanen = 'dl1(mss,a,deltas,dim_x)';  %pseudo-linear DL
%i.e. penalise for description of wieghts only
%Rissanen2= 'dl2(mss,unscaled_a,unscaled_deltas,dim_x,base,basic,offset,v)'; %approx DL
Rissanen2= 'dl2(mss,a,deltas,dim_x,base,basic,offset,v)'; %approx DL
%assume penalty for all significant parameters are equal and equal to the penalty for
%the weights (this assumption is not necessarily well founded --- just convenient)



%check for inputs
na=nargin;
if na<13,
  meth=[];
end;
if na<12,
   wobble=[];
end;
if na<11,
   nexpand=[];
end;
if na<10,
   nbasis=[];
end;
if na<9,
   nbuilds=[];
end;
if na<8,
   niterat=[];
end;
if na<7,
   penalty=[];
end;
if na<6,
   func=[];
end;
if na<5,
   v=[];
end;
if na<4,
   y_fit=[];
end;
if na<3,
  x_fit=[];
end;
if na<2,
   y_test=[];
end;
if na<1,
  x_test=[];
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
end;
if isempty(v),
   v=[-1 0 1 2];
   disp('v=');
   disp(v);
end;

%embed data, sort out y_test and y_fit
[dx,nx]=size(x_fit);
[dy,ny]=size(y_fit);
if dy~=1,
  disp('y_fit must be a row vector');
end;
if nx~=ny,
  disp('dimension of x_fit and length of y_fit don''t agree');
end;
uv=v(~isnan(v(:)));
uv=max(uv)+1;
if dx~=uv,
  disp('dimension of x_fit and v don''t agree');
end;
y=y_fit;
X=x_fit;
[dx,nx]=size(x_test);
[dy,ny]=size(y_test);
if dy~=1,
  disp('y_test must be a row vector');
end;
if nx~=ny,
  disp('dimension of x_test and length of y_test don''t agree');
end;
if dx~=uv,
  disp('dimension of x_test and v don''t agree');
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
      disp(['* Current best model : ',int2str(length(best_lambda(:,1))),' basis functions, RMS =',num2str(rms(best_err')),', MDL =',num2str(best_mdl),'.']);
      end;
   end;
end;


%build the model
for i=1:niterat
   base=pbase;
   for j=1:nbuilds,
      if trace,
         disp(['Starting calculation for node (',int2str(i),',',int2str(j),')']);
      end;
      [base,lambda,delta,basis,err,mdl]=rb_topdown(X,y_fit,v,func,base,penalty,nbasis,nexpand,wobble,meth);
      disp(['Calculation complete for node (',int2str(i),',',int2str(j),')']);
      %is latest model the best?
      if mdl<best_mdl
         improved=1;
         disp(['* Current best model : ',int2str(length(lambda(:,1))),' basis functions, RMS =',num2str(rms(err')),', MDL =',num2str(mdl),'.']);
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
         rb_set_globals(X,y_fit,best_base,best_lambda,best_delta,best_basis,penalty,best_mdl,best_err,v,func,meth);
      elseif trace,
         close;
      end;
   end;
end;

if ~improved,
   disp('Didn''t improve initial model --- bummer.');
end;



%proceed to do some testing.
y_test=[];
if isempty(y_test),
   %do nothing.
else,
   %one step predictions
   figure;
   %dishonest predictions
   [yt,yp,ep]=predict(x_fit,y_fit);
   subplot(211);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['dishonest predictions RMS = ',num2str(rms(ep))]);
   %honest predictions
   [yt,yp,ep]=predict(x_test,y_test);
   subplot(212);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['honest predictions RMS = ',num2str(rms(ep))]);
   %and a free run
   figure
   bounds=10*std(y_test)+abs(mean(y_test));
   %no-noise
   [yt,yp,ep]=freerun(y_test,1,min(2000,length(y_test)),[],bounds);
   subplot(211);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title('long term prediction');
   %noise
   [yt,yp,ep]=freerun(y_test,1,min(2000,length(y_test)),rms(best_err'),bounds);
   subplot(212);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['noisy simulation, std = ',num2str(rms(best_err'))]);
end;
