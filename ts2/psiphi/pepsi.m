%compare the performance of cubic hidden lasyer
% (cva) or radial basis (hl) upper layer to iterated one-step predictions.
%and draw a picture.
%
% for comparison of hidden layer radial basis network as well see pepsicola
%in /matlab/hidden
%

predlength=100;
len=1500;
predvect=[1 2 3 5 10 15 20 30 40];%[1 2 3 6 12 24 36 48 72];%[1 2 5 10 20 30 40 60 80];%[1 2 3 6 9 12 18 24 36];

[yp,yd]=cva_makeit(y(len:length(y)),predlength,1:predlength);  
e_pl=RMS(yp'-yd);

clear e_poly;

for np=1:3,
  [z_pred,y_pred,e_pred,y_data,B,Bsc,C,Csc] = ...
      cva_analysis(y,len,predlength,predvect,np,np);
  e_poly(:,np)=RMS(e_pred);
end;

close all
figure;axis([0 predlength 0 8000]);
plot(e_pl,'k');hold on;
plot(e_poly,'b--');
for i=1:length(predvect);
  plot(predvect(i)*[1 1],[0 8000],'g:');
end;
xlabel('prediction step');
ylabel('RMS(prediction error)');