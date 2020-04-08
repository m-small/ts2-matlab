function bftestmodel(y_fit,y_test);

%function bftestmodel(y_fit,y_test);
% or testmodel(y_test);
%
%
% M. Small 
% Created: 15/10/99
% Updated: 15/10/99

if nargin==1,
   y_test=y_fit;
   y_fit=[];
elseif max(size(y_test))==1,
   yn=y_test;
   y_test=y_fit(yn:end);
   y_fit=y_fit(1:(yn-1));
end;

global rb_x
[dx,nx]=size(rb_x);
mu=rb_x(dx,:);


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
