function fstr=rb_string(basic,base,offset,func,newold,deleted);
%
% M. Small 
% Created: 3/4/98
% Updated: 20/8/98

fstr='';
for i=1:length(basic);
   if basic(i)>offset,
      fnow=func{base.func(basic(i)-offset)};
      if newold(i),
         fnow=upper(fnow);
      end;
      fstr=[fstr,fnow(1)];
   elseif basic(i)==1,
      fstr=[fstr,'0'];
   elseif basic(i)<=offset,
      fstr=[fstr,'1'];
   else 
      basic(i)='-';      
   end;
end;

if ~isempty(deleted),
   fstr=[fstr(1:(deleted-1)),'#',fstr(deleted:end)];
end;
