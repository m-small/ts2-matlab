function combine(filename1,filename2);

%function combine(file1,file2);
% function combine(files);
%
%combine multiple pl_timeseries type models into a single ensemble model.
%the resultant model is stored as global variables
%
%either combine the two models stored in file1.mat and file2.mat or combine
%all models with filenames matching the pattern files
%

rb_get_globals;

if nargin<2,
    filename2=dir(filename1);
    if ~isempty(filename2),
        files=struct2cell(filename2);
        files=files(1,:);
    else,
        warning('Nothing found to combine');
        files=[];
    end;
else,
    files={filename1,filename2};
end;

nf=length(files);

cb_base.centres=[];
cb_base.strategy=[];
cb_base.func=[];
cb_base.radii=[];
cb_basis=[];
cd_delta=[];
cb_lambda=[];

ref_x=[];



tol=1e-9;% tolerance for independence

nbasis=0;

for filen=1:nf,
    load(files{filen});
    disp(['loading : ',files{filen}]);
    
    if isempty(ref_x),%check compatibility,
        ref_embed=rb_embed;
        ref_functions=rb_functions;
        ref_method=rb_method;
        ref_penalty=rb_penalty;
        ref_x=rb_x;
        ref_y=rb_y;

        %some other initialisation:
        %compute the true offset
        [nx,dx]=size(rb_x);
        offset=0;
        if any(rb_method=='c'),
            offset=offset+1;
        end;
        if any(rb_method=='l'),
            offset=offset+nx;
        end;
        cb_basis=1:offset;
        cb_lambda(cb_basis,1:2)=0;
        cb_delta(cb_basis,1)=0;
    else,
        if ~isequalwithequalnans(rb_embed,ref_embed),
            error('Embeddings not identical');
        end;
        if ~isequal(ref_functions,rb_functions),
            error('Candidate basis function not the same');
        end;
        if ~isequal(ref_method,rb_method),
            error('Model methodology not the same');
        end;
        if ~isequal(ref_penalty,rb_penalty),
            if ~isnan(ref_penalty),
                warning('Modelling penalty not uniform');
                ref_penalty=nan;
            end;
        end;
        if ~isequal(ref_x,rb_x) | ~isequal(ref_y,rb_y),
            warning('Data non-identical');
            ref_x=nan;
            ref_y=nan;
        end;
    end;



    %now concatenate
    nbasis=length(cb_base.strategy);
    cb_base.centres=[cb_base.centres rb_base.centres];
    cb_base.strategy=[cb_base.strategy; rb_base.strategy];
    cb_base.func=[cb_base.func; rb_base.func];
    cb_base.radii=[cb_base.radii rb_base.radii];
    cb_basis=[cb_basis rb_basis(rb_basis>offset)+nbasis];
    linearbits=find(rb_basis<=offset);
    cb_lambda(rb_basis(linearbits),:)=cb_lambda(rb_basis(linearbits),:)+rb_lambda(linearbits,:);
    cb_lambda=[cb_lambda;rb_lambda(rb_basis>offset,:)];
    cb_delta(rb_basis(linearbits))=cb_delta(rb_basis(linearbits))+rb_delta(linearbits);
    cb_delta=[cb_delta;rb_delta(rb_basis>offset)];
    
end;

%average
dontneed=find(cb_lambda(:,1)==0);
cb_lambda(dontneed,:)=[];
cb_delta(dontneed)=[];
cb_basis(dontneed)=[];

if isscalar(ref_x) & isnan(ref_x),
    %this should be sufficient... but it is not

    cb_lambda(:,2)=cb_lambda(:,2).*nf;
    ref_error=nan;
else
    %so we need to do this:
    phi=rb_Phi(ref_x,cb_base,ref_embed,ref_functions,ref_method); %compute phi
    phi=phi(:,cb_basis);
    [phi,scale]=normalize(phi); %and normalise it
    %remove redundant (linearly dependent) basis functions
    ntot=min(size(phi));
    disp('decomposing ...');
    [q,r,e]=qr(phi);
    disp(['rank = ',int2str(rank(phi))]);
    nout=sum(abs(diag(r))<tol);
    if nout>0,
        for i=1:nout,
            ind(i)=find(e(:,ntot-nout+i));
        end;
        disp(['eliminating ',int2str(nout),' nearly redundant basis functions']);
        cb_base.centres(:,cb_basis(ind)-offset)=[];
        cb_base.strategy(cb_basis(ind)-offset)=[];
        cb_base.func(cb_basis(ind)-offset)=[];
        cb_base.radii(:,cb_basis(ind)-offset)=[];
        cb_delta(cb_basis(ind)-offset)=[];        
        %the following line adjusts the indices beyond the redundant
        %components ... go figure!
        cb_basis=cb_basis-sum([zeros(1,ntot); ones(nout,1)*cb_basis>cb_basis(ind)'*ones(1,ntot) & ones(nout,1)*cb_basis>offset]);
        cb_basis(ind)=[];
       % cb_basis=cb_basis-sum([zeros(1,ntot-nout); ones(nout,1)*cb_basis>(ind+offset-1)'*ones(1,ntot-nout) & ones(nout,1)*cb_basis>offset]);
        % cb_basis=cb_basis-sum([zeros(1,ntot-nout); ones(nout,1)*(1:ntot-nout)>(ind+offset-1)'*ones(1,ntot-nout) & ones(nout,1)*cb_basis>offset]);
        %redo phi
        phi=rb_Phi(ref_x,cb_base,ref_embed,ref_functions,ref_method); %compute phi
        phi=phi(:,cb_basis);
        [phi,scale]=normalize(phi); %and normalise it
        %[q,r,e]=qr(phi);
    end;

    %recompute
    disp(['computing ',int2str(ntot-nout),' weights']);
    a=phi\ref_y';
    % a=r\(q'*ref_y');
    ref_error= ref_y' - phi*a;
    cb_lambda = [a scale'];% to avoid numerical error store normalized values

end;

%store
rb_base=cb_base;
rb_basis=cb_basis;
rb_delta=cb_delta;
rb_descr_length=nan;
rb_embed=ref_embed;
rb_error=ref_error;
rb_functions=ref_functions;
rb_lambda=cb_lambda;
rb_method=ref_method;
rb_penalty=ref_penalty;
rb_x=ref_x;
rb_y=ref_y;
