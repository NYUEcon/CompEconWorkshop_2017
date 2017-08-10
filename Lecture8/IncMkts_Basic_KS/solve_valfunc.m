function [v,jac] = solve_valfunc(c,s,param,glob,xxx)
%SOLVE_VALFUNC_NOAGG Solves value function in no uncertainty case
%-------------------------------------------------
%   Solves for the value function coeffivient vectors (cC,cE)
%
%   INPUTS
%   - c         = current collocation coefficient matrix
%   - bKS       = Krusell Smith algorithm coefficients.
%   - s         = state space
%   - param     = 
%   - glob      =
%   - options   = 
%   OUTPUT
%   - v         = 
%   - jac       = Jacobian of the value functions  
%-------------------------------------------------

%% Unpack coefficient vector
cv = c(1:end/2);
cE = c(end/2+1:end);   

%% solve value function for any (asset,labour productivity) state
B                       = menufun('bounds',s,[],param,glob);
obj                     = @(kstar)valfunc(cE,s,kstar,param,glob);
kstar                   = goldenx(obj,B(:,1),B(:,2));
[vV, Phi_kp_e_A_Kp]     = valfunc(cE,s,kstar,param,glob);

%% Compute vE and jacobian if requested
vE  = [];

Kprime     = menufun('Kprime', s, [], param, glob);
Kprime     = max(min(Kprime,max(glob.Kgrid0)),min(glob.Kgrid0));  % Make sure stays on grid
Phi_K      = splibas(glob.Kgrid0,0,glob.spliorder(4),Kprime);
Phi_k_e_A_Kp    = dprod(Phi_K, glob.Phi_keA);

if (nargin<=5)  
    % Expected value function
    vE = glob.Emat*Phi_k_e_A_Kp*cv;
end

if (nargout==2)
    jac = [ glob.Phi,           -param.beta*Phi_kp_e_A_Kp;
          - glob.Emat*Phi_k_e_A_Kp,           glob.Phi];          
end

%% Packup output

% Value function coefficients
v.vV       = vV;
v.vE       = vE;

% Optimal policy functions (correspond to each state element in 's').
v.kstar    = kstar;   % optimal asset holdings
v.kdist    = kstar;   % distribution of assets (next period) across agents



