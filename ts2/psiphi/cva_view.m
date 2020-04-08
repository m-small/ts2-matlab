function [img,imgx,imgy]=cva_view(yp,ep,yt,pl,len,A,Asc,ANp,n_obs,n_dyn,n_it,img_sz)

% [img,imgx,imgy] = ... 
%               cva_view(yp,ep,yt,pl,len,A,Asc,ANp,n_obs,n_dyn,n_it,img_sz)
%
% super fancy, extra special pictures.
% 
% This does the same thing as cva_slomo, and a bit more. A pdf of the
% prediction values for "dynamic" and "obseverational" noise n_dyn and
% n_obs. Dynamic noise is noise added at each step of the nonlinear (pl_)
% model predictions (to get the little red dots). The "observational" noise
% is noise on the data vector values (i.e. thge embedded vector point from
% which the future predictions are made), and is therefore a combination of
% what would conventionally be labelled dyn. and obs. noise.
%
% press p to generate the pdf (instruction on the effect of key presses
% appear in matlab window and on the figure title bar).
%
% Unfortunately there are several input arguments, they must be specified in
% the correct order to omit an argument specify it as []. The arguments are:
%    yp,ep,yt,pl --- predictions, errorbar predictions, test data and
% nonlinear (pl_) model predictions, all calculated by calls to the various
% functions in the cva directory.
%    len --- length of data to preview (i.e. upper bound on datum no. to
% calculate predictions from)
%    A,Asc --- A matrix (B matrix, whatever) the linear interpolation matrix
% and scaling factors.
%    ANp --- Order of polynomial interpolation/extrapolation
%    n_obs,n_dyn --- dynamic and observational noise levels, as descibed
% above.
%    n_it --- number of future predictions to estimate pdf from
%    img_sz --- number of bins (vertical) in histogram from which image of
% pdf is calculated.
%
% Outputs are:
%    img,imgx,imgy --- the last image calculated and the x and y ranges (NB:
% coloursc(imgx,imgy,img) will recreate this image, but because matlab draw
% images with the origin in the top left hand corner the YDir property of
% the current figure must be set to 'Normal').
%
% MAS 10-9-98

% Copyright (c) 1998 by Michael Small.
%
% Please see the copyright notice included in this distribution
% for full details.
%
%
% File   cva_movie.m
%   $Id$
%
% Created by Michael Small (<watchman@>) on Wed Aug 12 1998
%
% $Log$

img=[];imgx=[];imgy=[];
colourmap=colormap;lc=length(colourmap(:,1)); %lc should be 64...

ver=version;
ver=eval(ver(1));
if ver<5,
  black='c';
else
  black='k';
end;

na=nargin;
if na<4
  len=pl;
  pl=[];
  npl=0;
else,
  npl=1;
end;
if na<12,
  img_sz=50;
end;
if na<11,
  n_it=1000;
end;
if na<10,
  n_dyn=1;
end;
if na<10,
  n_obs=1;
end;

instruct='x - exit; p/q - calculate pdf (of cva/pl_model); b - back; other key - forward.';
figure(gcf);
set(gcf,'Name',instruct);disp(instruct);
clf;

miny=min([min(yt(:)),min(pl)]);
maxy=max([max(yt(:)),max(pl)]);
ps=length(yt(:,1));

global pl_pred_vect
ppv=pl_pred_vect;
ppv=-ppv;%(ppv<=0);

global pl_embed
es=max(pl_embed(~isnan(pl_embed) & ~ isinf(pl_embed)));

dv=0.95*es/(es+ps)+0.025;
predwin=subplot('position',[dv+0.05 0.11 0.925-dv 0.815]);
hold on;
datawin=subplot('position',[0.025 0.11 dv 0.815]);

axes(predwin);
%axis([0 ps+1 miny maxy]);
set(predwin,'DrawMode','fast');

axes(datawin);
%axis([-es-1 1 miny maxy]);
set(datawin,'DrawMode','fast');

%setup the initial picture
i=1;
axes(predwin);
pt=plot(yt(:,i),[black,'-.']);hold on;pp=plot(yp(:,i));
pv=plot(ppv+1,pl(i,:),'r*');  
set(pt,'EraseMode','background');
set(pp,'EraseMode','background');
set(pv,'EraseMode','background');
if npl,
  pph=plot(yp(:,i)+ep(:,i),'g-');
  ppl=plot(yp(:,i)-ep(:,i),'g-');
  set(pph,'EraseMode','background'); 
  set(ppl,'EraseMode','background');
  xlbl='data (black dot-dashed); cva-prediction (blue solid) and Std. Error (green solid); radial basis prediction (red stars)';
else
  xlbl='data (black dot-dashed); cva-prediction (blue solid); radial basis prediction (red stars)';
end;

hold off;title(['Test data (honest predictions) --- Datum:',int2str(i)]);
xlabel(xlbl);
drawnow;
axis([0 ps+1 miny maxy]);
axes(datawin);
pd=plot((-es):0,yt(1:(es+1),i),black);
set(pd,'EraseMode','background');
axis([-es-1 1 miny maxy]);
xlabel('data vector');hold off;
imgy=[maxy miny];
imgx=[2 ps];

%and up up, and away,
while i<=len,
  set(pt,'YData',yt(:,i));
  set(pp,'YData',yp(:,i));
  if npl,
    set(pph,'YData',yp(:,i)+ep(:,i));
    set(ppl,'YData',yp(:,i)-ep(:,i));
  end;
  set(pv,'YData',pl(i,:));
  set(get(predwin,'Title'),'String',['Test data (honest predictions) --- Datum:',int2str(i)]);
  drawnow;
  if i>es,
    set(pd,'YData',yt(1:(es+1),i-es));
  end;
  waitforbuttonpress;
  keypress=lower(get(gcf,'CurrentCharacter'));

  if keypress=='b'
    i=i-1;
    if i==0,
  	  i=1;
    end;
  
  elseif keypress=='x'
  
    break;
  
  elseif keypress=='p' | keypress=='q'

    if na<10,
      disp('Insufficient input args.');
    elseif i>es,
      ytp=flipud(yt(1:(es+1),i-es));
      if keypress=='q',
	[ypees,img]=cva_reddot(ytp,n_obs,n_dyn,n_it,ps,img_sz,[miny maxy]);
	titl='distribution of pl\_model predictions';
      else
	[ypees,img]=cva_pdf(ytp,n_obs,n_dyn,n_it,ppv,A,Asc,ANp,img_sz,[miny maxy]);
	titl='distribution of cva predictions';
      end;
      clear ypees;
      set(datawin,'color',colourmap(lc,:));
      axes(predwin);
      if ver>4,
        imghan=imagesc(imgx',imgy',flipud(img));
      else
	imghan=imagesc(imgx',imgy',img);
      end;
      title([titl,' --- Datum:',int2str(i)]);
      xlabel(xlbl);
      set(get(imghan,'Parent'),'YDir','normal'); %imagesc is arse over.
      hold on;
      pt=plot(yt(:,i),[black,'-']);pp=plot(yp(:,i),'r');
      pv=plot(ppv+1,pl(i,:),'r*');  
      if npl,
        pph=plot(yp(:,i)+ep(:,i),'y-');
        ppl=plot(yp(:,i)-ep(:,i),'y-');hold off;
      end;
      waitforbuttonpress;
      if lower(get(gcf,'CurrentCharacter'))=='x'
	break
      end;
      cla;
      set(predwin,'DrawMode','fast'); hold on;
      title(['Test data (honest predictions) --- Datum:',int2str(i)]);
      xlabel(xlbl);
      pt=plot(yt(:,i),[black,'-.']);pp=plot(yp(:,i));
      pv=plot(ppv+1,pl(i,:),'r*');  
      set(pt,'EraseMode','background');
      set(pp,'EraseMode','background');
      set(pv,'EraseMode','background');
      if npl,
        pph=plot(yp(:,i)+ep(:,i),'g-');
        ppl=plot(yp(:,i)-ep(:,i),'g-');hold off;      
        set(pph,'EraseMode','background'); 
        set(ppl,'EraseMode','background');
      end;
      
      axes(datawin);cla;
      set(datawin,'DrawMode','fast');
      pd=plot((-es):0,yt(1:(es+1),i),black);
      set(pd,'EraseMode','background');
      axis([-es-1 1 miny maxy]);
      xlabel('data vector');hold off;

    else,
      disp(['Insufficient data, wait for datum ',int2str(es+1)]);
    end;

  else,
    i=i+1;
  end;     
end;

% End of cva_view.m





