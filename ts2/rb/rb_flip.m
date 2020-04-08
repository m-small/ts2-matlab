function rb_flip(basename,inout)

%function rb_flip(basename,inout)
%
%move the rb_model storted in global variables [basename,'_rb_*'] to 'rb_*'
% (if inout=1) or viceversa (inout=0)

rb_get_globals
varname={'rb_x','rb_y','rb_base','rb_lambda','rb_delta','rb_basis','rb_penalty','rb_descr_length','rb_error','rb_embed','rb_functions','rb_timescale'};
nv=length(varname);

for i=1:nv,
    eval(['global ',basename,'_',varname{i}]);
end;

switch inout,
    case 1
        for i=1:nv,
            eval([varname{i},'=',basename,'_',varname{i},';']);
        end;
    case 0
        for i=1:nv,
            eval([basename,'_',varname{i},'=',varname{i},';']);
        end;
    case -1
        for i=1:nv,
            eval(['clear global ',basename,'_',varname{i},';']);
        end;
 
    otherwise
        disp('buggered');
end;
