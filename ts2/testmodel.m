function testmodel(y_fit,y_test);

%function testmodel(y_fit,y_test);
% or testmodel(y_test);
%
%
% M. Small 
% Created: 15/6/99
% Updated: 15/6/99

if nargin==1,
   y_test=y_fit;
   y_fit=[];
elseif max(size(y_test))==1,
   yn=y_test;
   y_test=y_fit(yn:end);
   y_fit=y_fit(1:(yn-1));
end;




%proceed to do some testing.
if isempty(y_test),
   %do nothing.
else,
   %one step predictions
   figure;
   %dishonest predictions
   [yt,yp,ep]=predict(y_fit);
   subplot(211);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   title(['dishonest predictions RMS = ',num2str(rms(ep))]);
   %honest predictions
   [yt,yp,ep]=predict(y_test);
   subplot(212);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on;
   errs=rms(ep);
   title(['honest predictions RMS = ',num2str(errs)]);
      xlabel('true (black), predicted (red), and error (blue)');
%and a free run
   figure
   bounds=1000*std(y_test)+abs(mean(y_test));
   %no-noise
   [yt,yp,ep]=freerun(y_test,1,min(2000,length(y_test)),0,bounds);
   subplot(211);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
   zoom on; drawnow;
   %noise
   [yt,yp,ep]=freerun(y_test,1,min(2000,length(y_test)),errs,bounds);
   subplot(212);
   plotcols(ep,'b:',yt,'k-',yp,'r--');
    xlabel('true (black), predicted (red), and error (blue)');
  zoom on;
end;
