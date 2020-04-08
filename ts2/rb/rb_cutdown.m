function rb_cutdown(x,y,insampl);

%function rb_cutdown(x,y,insampl);
% 
%Cut a rb_model down to size ala cross validation
% rb_eval(x)~y
%the first insampl observation are fit data, the resst is test data
%

rb_get_globals;

% extract model data
[dc,nc]= size(rb_base.centres); 				% centres tell us dimension

%assume x and y are compatible and insampl<length(y)
[dx,nx]= size(x);
y=y(:);
ly=length(y);

v= rb_embed;
[nv,dv]= size(v);

X=rb_Phi(x,rb_base,v,rb_functions,rb_method);

ll=length(rb_basis);
[dX,lX]=size(X);
if ll~=length(rb_basis), 
 disp(['WARNING: ll~=lX, ll=',int2str(ll),'; lX=',int2str(lX),'.']);
 %X= X(:,1:ll);
end;

fiterr=[];
testerr=[];
for nb=1:length(rb_basis),
   
   lambdas=X(1:insampl,rb_basis(1:nb))\y(1:insampl);
   yp= X(:,rb_basis(1:nb)) * lambdas;
   
   fiterr(nb)=rms(yp(1:insampl)'-y(1:insampl)');
   testerr(nb)=rms(yp((insampl+1):end)'-y((insampl+1):end)');
   
   if trace,
      subplot(311);
      plot(yp(1:insampl),'b');
      hold on;
      plot(y(1:insampl),'y');
      hold off;
      title(['In-sample (nb=',int2str(nb),' of ',int2str(length(rb_basis)),')']);
      subplot(312);
      plot(yp((insampl+1):end),'b');
      hold on;
      plot(y((insampl+1):end),'y');
      hold off;
      title('Out-sample (data (blue), prediction (y))');
      subplot(313);
      plot(fiterr,'g');
      hold on;
      plot(testerr,'c');
      hold off;
      title('RMS error: fit (green), test (blue)');
      xlabel('model size');
      drawnow;
   end;
end;

[merr,msize]=min(testerr);
if trace,
   subplot(313);
   hold on;
   axis('tight');
   ax=axis;
   plot(msize,merr,'ro');
   plot(msize*[1 1],ax(1:2)','r-.');
   hold off;
end;

rb_basis=rb_basis(1:msize);
rb_delta=rb_delta(1:msize);
rb_lambda=rb_lambda(1:msize,:);

