function [Pz,zvec,Pssz] = setup_Markov(Nz,sige,rhoz,A)
%________________________________________________________________________
% Function
% 1. Using rouwenhorst
sig_uncond      = sige/sqrt(1-rhoz^2);
mu_uncond       = log(A)-(1/2)*(sige^2)/((1+rhoz)*(1-rhoz));
[Pz, logzvec]   = setup_rouwen(rhoz,mu_uncond,sig_uncond,Nz);
Pz              = Pz';
zvec            = exp(logzvec);
Pssz            = Pz^10000;
Pssz            = Pssz(1,:)';
% 2. Remove points with Pr(z)<1e-5 in the stationary distribution
% Issue: The rouwen procedure creates many points with little prob
cut = 1e-6;       % Chop off points with low productivity
if (cut>1e-8)&&(min(Pssz)<cut)
    go  = 1;
    Nztemp = Nz;
    while go
        Nztemp          = Nztemp+1;
        [Pz, logzvec]   = setup_rouwen(rhoz,mu_uncond,sig_uncond,Nztemp);
        Pz              = Pz';
        zvec            = exp(logzvec);
        Pssz            = Pz^10000;
        Pssz            = Pssz(1,:)';
        X               = sum(Pssz>cut);
        if X>=Nz;go=0;end;
    end
    i       = (Pssz>cut);
    zvec    = zvec(i);
    Pz      = Pz(i,i);
end
Pz          = bsxfun(@rdivide,Pz,sum(Pz,2));
Pz(:,end)   = 1-sum(Pz(:,1:end-1),2);
Pssz        = Pz^10000;
Pssz        = Pssz(1,:)';

