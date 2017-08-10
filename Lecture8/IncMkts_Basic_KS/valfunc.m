function [ v, Phi_keAK ] = valfunc(cE,s,kprime,param,glob)
%VALFUNC gives value function value given parameters, states
%-------------------------------------------------
%   Computes the the value function
%
%   INPUTS
%   - cE        = Collocation coefficient matrix for expected value function
%   - s         = State space
%   - kprime    = Optimal individual capital choice
%   - param     = Model parameters
%   - glob      = Model global variables
%   OUTPUT
%   - v         = Value function 
%   - Phi_keAK  = Basis matrix updated for capital choice, kprime
    
%-------------------------------------------------

%% Compute felicity function + continuation value

F           = menufun('felicity',s,kprime,param,glob);
Phi_k       = splibas(glob.kgrid0,0,glob.spliorder(1),kprime);

% Kprime     = menufun('Kprime', s, [], param, glob);
% Kprime     = max(min(Kprime,max(glob.Kgrid)),min(glob.Kgrid));  % Make sure stays on grid

% Phi_K      = splibas(glob.Kgrid0,0,glob.spliorder(4),Kprime);
Phi_keAK    = dprod(glob. Phi_K, dprod(glob.Phi_eA, Phi_k));

% Phi_keAK    = dprod(glob.Phi_eAK, Phi_k);

v           = F + param.beta*Phi_keAK*cE;


end
