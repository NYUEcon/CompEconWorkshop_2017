function [ eq ] = solve_statdist(ck,K,param,glob,options)
%SOLVE_STATDIST solves for the stationary distribution of the model
%-------------------------------------------------
%
%   INPUTS
%   - c         = collocation coefficients
%   - K         = current guess for aggregate capital
%   - param     = model parameters
%   - glob      = global variables
%   - options   = plotting options
%
%   OUTPUT
%   - eq        = outputs for the stationary distribution
%-------------------------------------------------

%% Compute stationary distribution
% Solve again on a finer grid for pP
kgrid           = glob.kgrid;
% sf              = glob.sf;
% glob.Phi_E      = glob.Phi_Ef;
% glob.Phi        = glob.Phif;
kdist         = glob.Phif*ck;
% v               = solve_valfunc(ck,sf,K,param,glob,1); 

% Compute stationary distribution
kdist           = min(kdist,max(kgrid));
fspaceergk      = fundef({'spli',glob.kgridf,0,1});
Qk              = funbas(fspaceergk,kdist);
QE              = glob.QE;
Q               = dprod(QE,Qk);

L               = ones(size(Q,1),1);
L               = L/sum(L);

for itL = (1:options.itermaxL);
    Lnew    = Q'*L;
    dL      = max(abs(Lnew-L));
    if (dL<options.tolL)
        break
    end
    
    if mod(itL,100)==0
        if strcmp(options.print,'Y')
            fprintf('dL:\t%1.3e\n',dL);
        end
    end
    L       = Lnew;
end


%% Compute aggregates and implied interest rate and wages

Kagg = L'*kdist;

%% Pack-up output
eq.L    = L;
eq.Q    = Q;
eq.Kagg = Kagg;



end