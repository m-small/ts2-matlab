function z=powerlaw(m,n,gamma)

%function z=powerlaw(m,n,gamma)
%function z=powerlaw(m,gamma)
%function z=powerlaw(gamma)
%
%m-by-n (or m-by-1, or 1-by-1 random sample of the power-law distribution
%with exponent gamma - p(x)=kx^-gamma (where k is 1/zeta(gamma-1) - or
%something)
%
%z=gprnd(1/gamma,1/gamma,1)
%for zeta ditribution need to round.
%


na=nargin;

switch na
    case 1
        %only one argument - it's gamma
        z=powerlaw(1,1,m);
    case 2
        %two args. - it's the vector length ... and gamma
        z=powerlaw(m,1,n);
    case 3,
        %we've got three arguments - m,n and gamma
        %disp(['gamma = ',num2str(gamma)]);
        gamma=abs(gamma); %otherwise it's divergent
        %pareto distribution with shape parameter 1/gamma and location x0=1
        z=gprnd(1/gamma,1/gamma,1*ones(m,n));
    otherwise
        z=powerlaw(2);
end;


