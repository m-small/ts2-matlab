function bifurcation(y_fit,y_test,v,func,penalty,niterat,nbuilds,nbasis,nexpand,wobble,meth);

%function bifurcation(y_fit,y_test,v,func,penalty,niterat,nbuilds,nbasis,nexpand,wobble,meth);
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
%  meth='clr'
%
% Build a radial basis model of a embedded time series, and include a bifurcation parameter 
% (bifurcation param. varies monotonically from minimum to maximum). v should include reference 
% to the bifurcation co-ord as infs.
%
% M. Small 
% Created: 6/6/99
% Updated: 15/6/00


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
   v=[];
end;
if na<2,
   y_fit=[];
end;
if na<1,
   y_test=[];
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
uv=v(~isnan(v(:)) & ~isinf(v(:)));
uv=min(uv):max(uv);
[X,y]=embed(y_fit,uv);
if length(y_test)==1,
   if round(y_test)==y_test,
      yy=y_test;
      X=X(:,1:(yy-1));
      y_test=y(yy:end);
      y=y(1:(yy-1));
      clear yy;      
   end;
end;
y_fit=y;

%now add the monotonic bifurcation param.and correct v
if isempty(X),
  mu=1:length(y);
  mu=(mu-1)/(length(mu)-1);
  X=mu;
else,
  mu=min(X(:)):((max(X(:))-min(X(:)))/(length(y)-1)):max(X(:));
  X=[X;mu];
end;
v(isinf(v))=max(uv)+1;
clear y

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
if isempty(y_test),
   %do nothing.
else,
   %one step predictions
   figure;
   %dishonest predictions
   [yt,yp,ep]=bfpredict(y_fit,mu);
   subplot(211);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['dishonest predictions RMS = ',num2str(rms(ep))]);
   %honest predictions
   mutest=(mu(end)+diff(mu(1:2))):diff(mu(1:2)):(mu(end)+diff(mu(1:2))*length(y_test));
   [yt,yp,ep]=bfpredict(y_test,mutest);
   subplot(212);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['honest predictions (extrapolation of bifurcation parameter) RMS = ',num2str(rms(ep))]);
   %and a free run
   figure
   bounds=10*std(y_test)+abs(mean(y_test));
   %no-noise
   mutest=[mu(1) mu(2)-mu(1)];
   [yt,yp,ep]=bffreerun(y_test,1,min(2000,length(y_test)),[],bounds,mutest);
   subplot(211);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title('long term prediction');
   %noise
   [yt,yp,ep]=bffreerun(y_test,1,min(2000,length(y_test)),rms(best_err'),bounds,mutest);
   subplot(212);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['noisy simulation, std = ',num2str(rms(best_err'))]);
end;
