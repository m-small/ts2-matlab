function [e_fit,e_test,ep]=testio(x_fit,y_fit,td);

%function testio(x,y,testdatum);
% or testio(x,y);
%
%
% M. Small 
% Created: 15/6/99
% Updated: 15/6/99

%embed data, sort out y_test and y_fit
[dx,nx]=size(x_fit);
y_fit=y_fit(:)';
[dy,ny]=size(y_fit);
if ny~=nx,
   error('size of x_fit and y_fit don''t match');
end;
if isempty(td),
   td=floor(nx/2);
end;
if td>nx,
   error('testdatum too big');
end;
y_test=y_fit((td+1):end);
x_test=x_fit(:,(td+1):end);
y_fit=y_fit(1:td);
x_fit=x_fit(:,1:td);

 %one step predictions
   figure;
   %dishonest predictions
yp=rb_eval(x_fit)';
ep=yp-y_fit;
yt=y_fit;   
subplot(211);
plotcols(ep,'b:',yt,'y-',yp,'r--');
zoom on;
title(['dishonest predictions RMS = ',num2str(rms(ep))]);

e_fit=rms(ep);

%honest predictions

yp=rb_eval(x_test)';
ep=yp-y_test;
yt=y_test;   
subplot(212);
plotcols(ep,'b:',yt,'y-',yp,'r--');
zoom on;
title(['honest predictions RMS = ',num2str(rms(ep))]);
xlabel('true (yellow), predicted (red), and error (blue)');


e_test=rms(ep);