function lag=findlag(y),

%make a "good" and robust (?) guess at an appropriate embedding lag

y=y(:);
ly=length(y);

n=min(max(100,ly/20),ly/2);

lag=inf;
while n<ly,
    [r,t]=acorr(y,n);
    if all(r>0),
        n=n*2;
        %all positive, no zeros
    else,
        r=r((length(r)+1)/2:end);
        [~,lag]=find(r<0);
        lag=min(lag)-1;
        n=inf;
    end;    
end;

if isinf(lag) %no zeros - first minimum?
    dd=diff(lag),
    [,lag]=min(dd>2);
    if ~isempty(lag),
        lag=lag(1)-1;
    else,
        lag=inf;
    end;
end;

if isinf(lag) % no zeros, no minimum - lag=1?
    lag=1;
end;

    