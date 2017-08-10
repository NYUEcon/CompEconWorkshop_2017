function [irfsim] = irf(prevdate, enddate, param,glob,options)
%IRF Plot impulse response functions 
%-------------------------------------------------
%
%   INPUTS
%   - param     = model parameters
%   - glob      = global variables
%   - options   = plotting options
%
%   OUTPUT
%   - irfsim    = Contains the IRF series

%-------------------------------------------------


% % Aggregate state
A_t = glob.Agrid(1)*ones(options.T,1);  % Begin in bad state
A_t(options.burn:options.burn+10) = glob.Agrid(2);   % Good state for 10 periods

% Initialize capital stock
K_t = zeros(options.T,1);
K_t(1,1) = param.Kss;   % Start aggregate capital at RA steady state

K_plm_t = zeros(options.T,1);
K_plm_t(1,1) = param.Kss;

% Initialize distribution
L_tmp = zeros(glob.Nkf, glob.Ne);
ix = (glob.kmax - glob.kmin)/glob.Nkf;
ix = round(param.Kss/ix);
L_tmp(ix, 1:2) = 1/2;

L_t = zeros(glob.Nkf*glob.Ne, options.T);
L_t(:,1) = reshape(L_tmp, glob.Nkf*glob.Ne,1);

kdist_t = zeros(glob.Nkf*glob.Ne, options.T);
cdist_t = zeros(glob.Nkf*glob.Ne, options.T);


sig_k_t = zeros(options.T,1);
% sig_ky_t = zeros(options.T,1);

Y_t = zeros(options.T,1);
C_t = zeros(options.T,1);
I_t = zeros(options.T,1);
r_t = zeros(options.T,1);
w_t = zeros(options.T,1);

% Solve individual problem for current KS parameters
solved = solve_coeff(param,glob,options);

for tt=1:options.T
    
    sf_t = gridmake(glob.kgridf, glob.egridf, A_t(tt), K_t(tt));
    
    Phi_k      = splibas(glob.kgrid0,0,glob.spliorder(1),sf_t(:,1));
    Phi_e      = splibas(glob.egrid0,0,glob.spliorder(2),sf_t(:,2));
    Phi_A      = splibas(glob.Agrid0,0,glob.spliorder(3),sf_t(:,3));
    Phi_K      = splibas(glob.Kgrid0,0,glob.spliorder(4),sf_t(:,4));
    Phi_keAK   = dprod(Phi_K, dprod(Phi_A, dprod(Phi_e, Phi_k)));
    
    % Get non-state aggregate variables
    Y_t(tt,:) = menufun('output',sf_t(1,:),[],param,glob);
    r_t(tt,:) = menufun('interest',sf_t(1,:),[],param,glob);
    w_t(tt,:) = menufun('wage',sf_t(1,:),[],param,glob);
    
    kdist_t(:,tt) = Phi_keAK*solved.ck; 
    kdist_t(:,tt)    = max( min(kdist_t(:,tt), max(glob.kgridf)), min(glob.kgridf));
    cdist_t(:,tt) = Phi_keAK*solved.cc; 
    
    
    % Compute stationary distribution
    fspaceergk      = fundef({'spli',glob.kgridf,0,1});
    Qk              = funbas(fspaceergk, kdist_t(:,tt));
    P_e             = (1-A_t(tt))*glob.Pe_b + A_t(tt)*glob.Pe_g;
    Qe              = kron(P_e, ones(glob.Nkf,1));
    Q_t             = dprod(Qe,Qk);
    % ^^^ NOTE: need to use the correct probability matrix for
    % the idiosyncratic state. It depends on the current aggregate
    % state, as that affects how likely the employment transitions are. 
    

    L_t(:,tt+1)     = Q_t'*L_t(:,tt);

    % NOTE: tomorrow's capitals depends on the CURRENT distribution over agents
    % making the choice of k_t+1. 
    K_t(tt+1) = L_t(:,tt)'*kdist_t(:,tt);
    K_t(tt+1) = max( min(K_t(tt+1), max(glob.Kgrid)), min(glob.Kgrid));
    
    C_t(tt) = L_t(:,tt)'*cdist_t(:,tt);
    I_t(tt) = Y_t(tt,:) - C_t(tt);
    
    sig_k_t(tt+1) = std(kdist_t(:,tt), L_t(:,tt));

    % Capital path according to PLM
    K_plm_t(tt+1) = exp( (1-A_t(tt))*(param.b_KS(1) + param.b_KS(2).*log(K_t(tt))) ...
        +  A_t(tt)*(param.b_KS(3) + param.b_KS(4).*log(K_t(tt))) );

   
end

% Plot IRFS
figure
subplot(2,3,1)
prod_t = (1-A_t)*(1-param.Da) + A_t*(1+param.Da);
mu = mean(prod_t(options.burn-prevdate:options.burn-1));
plot((prod_t(options.burn-prevdate:options.burn+enddate)-mu)/mu,'linewidth',2)
title('TFP')
grid on

subplot(2,3,2)
mu = mean(Y_t(options.burn-prevdate:options.burn-1));
plot((Y_t(options.burn-prevdate:options.burn+enddate) - mu)/mu,'linewidth',2)
grid on
title('Aggregate output')

subplot(2,3,3)
mu = mean(K_t(options.burn-prevdate:options.burn-1));
plot((K_t(options.burn-prevdate:options.burn+enddate) - mu)/mu,'linewidth',2)
grid on
title('Aggregate capital stock')

subplot(2,3,4)
mu = mean(C_t(options.burn-prevdate:options.burn-1));
plot((C_t(options.burn-prevdate:options.burn+enddate) - mu)/mu,'linewidth',2)
grid on
title('Aggregate Consumption')

subplot(2,3,5)
mu = mean(I_t(options.burn-prevdate:options.burn-1));
plot((I_t(options.burn-prevdate:options.burn+enddate) - mu)/mu,'linewidth',2)
grid on
title('Aggregate investment')

subplot(2,3,6)
mu = mean(sig_k_t(options.burn-prevdate:options.burn-1));
plot((sig_k_t(options.burn-prevdate:options.burn+enddate) - mu)/mu,'linewidth',2)
title('Std. dev. of wealth')
grid on



% Pack up output:
irfsim.A_t = A_t;
irfsim.K_t = K_t; 
irfsim.K_plm_t = K_plm_t;
irfsim.Y_t = Y_t;
irfsim.C_t = C_t;
irfsim.I_t = I_t;
irfsim.r_t = r_t; 
irfsim.w_t = w_t; 
irfsim.L_t = L_t;
irfsim.kdist_t = kdist_t;
irfsim.sig_k_t = sig_k_t;






end