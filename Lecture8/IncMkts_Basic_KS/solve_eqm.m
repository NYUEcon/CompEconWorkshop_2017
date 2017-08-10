function [ eq ] = solve_eqm(param,glob,options)
%SOLVE_EQM solves for the equilibrium interest rate
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

% Initial guess for aggregate capital stock
K = param.kss;

damp = 1;
for iteq = (1:options.itermaxKeq);
    
    solved = solve_coeff(K, param,glob,options);
    ck = solved.ck;
    eq = solve_statdist(ck,K,param,glob,options);
    dK  = abs(eq.Kagg-K)/K;
    if (dK<options.tolKeq)
        break
    end
    
    if strcmp(options.eqmprint,'Y')
        fprintf('Iteration = \t%i: K=\t%1.4f  A=\t%1.4f\n',iteq, K, eq.Kagg)
    end
    % Update aggregate capital using dampening
    K = damp*K + (1-damp)*eq.Kagg;
    damp = 0.9995*damp;
end

% Note the damping procedure: the damping parameter shifts more weight
% to updating as the algorithm iterates. This means the algorithm is slow
% to incoporate new information, however it prevents wild oscillations in 
% capital while trying to find the fixed point. This is particularly useful
% if the functions (capital demand and supply) are not Lipschitz bounded in
% a neighborhood of the initial condition. This is indeed the case, because
% when capital is below the RA steady state, the interest rate is greater
% than 1/beta - 1 + delta, which induces non-stationarity in asset demand 
% by agents. 


%% Compute equilibrium objects. 
Y = menufun('output',[],[],K,param,glob);
r = menufun('interest',[],[],K,param,glob);
w = menufun('wage',[],[],K,param,glob);

%% Pack-up output
eq.K = K;
eq.Y = Y;
eq.r = r;
eq.w = w;
eq.solved = solved;



end