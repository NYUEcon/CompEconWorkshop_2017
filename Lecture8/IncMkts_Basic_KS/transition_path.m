function [ eq ] = transition_path(ck,K,param,glob,options)
%TRANSITION_PATH solves for the transition path following exogenous shock
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
options.print = 'N';
options.eqmprint = 'Y';

T = 150;
damp = 0.9995;
K = param.kss*ones(T,1);    % Guess a sequence of K
A = zeros(T,1);             % Actual assets
s = glob.s;

% Simulate exogenous shock path to depreciation
rho = 0.75;
delta = zeros(T,1);
eps = zeros(T,1);
delta(1) = param.delta;
eps(2) = 0.025;
for tt = 2:T
   delta(tt) = (1-rho)*param.delta + rho*delta(tt-1) + eps(tt); 
end

% Solve the date t=1 problem
param.delta = delta(1);
eq = solve_eqm(param,glob,options);
K(1) = eq.K;
A(1) = eq.K;

% Solve the date t=T problem
param.delta = delta(T);
% eq = solve_eqm(param,glob,options);
K(T) = eq.K;
A(T) = eq.K;
solved = solve_coeff(K(T), param, glob,options);
cE = solved.c(end/2+1:end);   % yields cE(T)

figure
for iter = 1:100

    % Now recursively solve problem from t=T-1 to t=2
    for tt = 2:T-1
        
        param.delta = delta(T+1-tt);
        B = menufun('bounds',s,[],K(T+1-tt),param,glob);
        obj = @(kstar)valfunc(cE,s,kstar,K(T+1-tt),param,glob);
        kstar = goldenx(obj,B(:,1),B(:,2));
        vV = valfunc(cE,s,kstar,K(T+1-tt),param,glob);
        vE = glob.Emat*vV;
        cE = glob.Phi\vE;     % yields cE(T+1-tt),  used as previous period's expectations
        ck = glob.Phi\kstar;  % yields cK(T+1-tt)
        
        eq = solve_statdist(ck,K(T+1-tt),param,glob,options);
        A(T+1-tt) = eq.Kagg;     % yields actual asset holdings in T+1-tt
        
    end
    
    plot(K, 'linewidth', 2)
    hold all
    plot(A,'linewidth', 2)
    hold off
    grid on
    legend('Capital', 'Aggregated capital')
    
    drawnow
    
    % Update aggregate capital using dampening
    damp = 0.95;
    K(2:T-1) = damp*K(2:T-1) + (1-damp)*A(2:T-1);
%     damp = 0.9995*damp;
    fprintf('Iteration = %i\n', iter)
end






end