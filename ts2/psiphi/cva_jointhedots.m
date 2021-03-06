function cva_jointhedots(yp,ep,yt,pl,len)

% =cva_jointhedots(yp,ep,yt,pl,len)  draw some pretty pictures
%
% same as cva_slomo but it there is a dotted red line connecting the red dots.


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

ver=version;
if eval(ver(1))<5,
  black='c';
else
  black='k';
end;

na=nargin;
if na<4
  len=pl;
  pl=[];
end;

figure(gcf);
clf;

miny=min(yt(:));%min([min(yt(:)),min(pl)]);
maxy=max(yt(:));%max([max(yt(:)),max(pl)]);
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

if na==5,
  disp('any key - go forward; ''b'' - go backward; ''x'' exit.');
  if isempty(ep),
    i=1;
    axes(predwin);
    pt=plot(yt(:,i),[black,'-.']);hold on;pp=plot(yp(:,i));
    pv=plot(ppv+1,pl(i,ppv+1),'r*');  
    pv2=plot(pl(i,:),'r:');  
    set(pt,'EraseMode','background');
    set(pp,'EraseMode','background');
    set(pv,'EraseMode','background');
    set(pv2,'EraseMode','background');
    hold off;
    title(['Test data (honest predictions) --- Datum:',int2str(i)]);
    xlabel('data (black dot-dashed); cva-prediction (blue solid); radial basis prediction (red stars)');
    drawnow;axis([0 ps+1 miny maxy]);
    axes(datawin);
    pd=plot((-es):0,yt(1:(es+1),i),black);
    set(pd,'EraseMode','background');
    axis([-es-1 1 miny maxy]);
    xlabel('data vector');hold off;
  while i<=len,
    set(pt,'YData',yt(:,i));
    set(pp,'YData',yp(:,i));
    set(pv2,'YData',pl(i,:));
    set(pv,'YData',pl(i,ppv+1));
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
    else,
      i=i+1;
    end;     
  end;
  
  else

    i=1;
    axes(predwin);
    pt=plot(yt(:,i),[black,'-.']);hold on;pp=plot(yp(:,i));
    pph=plot(yp(:,i)+ep(:,i),'g-');
    ppl=plot(yp(:,i)-ep(:,i),'g-');
    pv=plot(ppv+1,pl(i,ppv+1),'r*');  
    pv2=plot(pl(i,:),'r:');  
    set(pt,'EraseMode','background');
    set(pp,'EraseMode','background');
    set(pph,'EraseMode','background');
    set(ppl,'EraseMode','background');
    set(pv,'EraseMode','background');
    set(pv2,'EraseMode','background');
    hold off;title(['Test data (honest predictions) --- Datum:',int2str(i)]);
    xlabel('data (black dot-dashed); cva-prediction (blue solid) and Std. Error (green solid); radial basis prediction (red stars)');
    drawnow;
    axis([0 ps+1 miny maxy]);
    axes(datawin);
    pd=plot((-es):0,yt(1:(es+1),i),black);
    set(pd,'EraseMode','background');
    axis([-es-1 1 miny maxy]);
    xlabel('data vector');hold off;
  while i<=len,
    set(pt,'YData',yt(:,i));
    set(pp,'YData',yp(:,i));
    set(pph,'YData',yp(:,i)+ep(:,i));
    set(ppl,'YData',yp(:,i)-ep(:,i));
    set(pv,'YData',pl(i,ppv+1));
    set(pv2,'YData',pl(i,:));
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
    else,
      i=i+1;
    end;     
  end;
end;

elseif na==4,

  i=1;
    axes(predwin);
    pt=plot(yt(:,i),[black,'-.']);hold on;pp=plot(yp(:,i));
    set(pt,'EraseMode','background');
    set(pp,'EraseMode','background');
    hold off;title(['Test data (honest predictions) --- Datum:',int2str(i)]);drawnow;
    axes(datawin);
    pd=plot((-es):0,yt(1:(es+1),i),black);
    set(pd,'EraseMode','background');
    axis([-es-1 1 miny maxy]);
    xlabel('data vector');hold off;
  for i=1:len,
    set(pt,'YData',yt(:,i));
    set(pp,'YData',yp(:,i));
    set(get(predwin,'Title'),'String',['Test data (honest predictions) --- Datum:',int2str(i)]);
    drawnow;
    if i>es,
      set(pd,'YData',yt(1:(es+1),i-es));
    end;
    waitforbuttonpress;
  end;
end;

% End of cva_movie.m





