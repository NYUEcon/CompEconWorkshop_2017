function [sim] = simulation(param,glob,options)
%SIMULATON Simulate the model
%-------------------------------------------------
%
%   INPUTS
%   - param     = model parameters
%   - glob      = global variables
%   - options   = plotting options
%
%   OUTPUT
%   - b_KS      = solved KS coefficients 
%   - R2bad     = R^2 for bad states
%   - R2good    = R^2 for good states
%-------------------------------------------------

% % Simulate aggregate shock
A_t = zeros(options.T,1);
A_t(1,1) = 0;     % Begin in bad state

rng(1);              % Make sure we are getting same simulation each time
tmp = rand(options.T,1);

for tt=2:options.T 
   if tmp(tt) <= glob.PA(1+A_t(tt-1), 1) 
      A_t(tt) = 0; % 
   else
      A_t(tt) = 1;
   end
end


% Initialize capital stock
K_t = zeros(options.T,1);
K_t(1,1) = param.Kss;   % Start aggregate capital at RA steady state

K_plm_t = zeros(options.T,1);
K_plm_t(1,1) = param.Kss;

K_plm0_t = zeros(options.T,1);
K_plm0_t(1,1) = param.Kss;

% Initialize distribution

L_tmp = zeros(glob.Nkf, glob.Ne);
ix = (glob.kmax - glob.kmin)/glob.Nkf;
ix = round(param.Kss/ix);
L_tmp(ix, 1:2) = 1/2;

L_t = zeros(glob.Nkf*glob.Ne, options.T);
L_t(:,1) = reshape(L_tmp, glob.Nkf*glob.Ne,1);



kdist_t = zeros(glob.Nkf*glob.Ne, options.T);

sig_k_t = zeros(options.T,1);
% sig_ky_t = zeros(options.T,1);

Y_t = zeros(options.T,1);
% C_t = zeros(options.T,1);
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
    
%     kdist_t(:,tt)    = min(kdist_t(:,tt),max(glob.kgridf));
    kdist_t(:,tt)    = max( min(kdist_t(:,tt), max(glob.kgridf)), min(glob.kgridf));
        
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
    K_t(tt+1) = L_t(:,tt)'*kdist_t(:,tt);    % L_t(:,tt) 
    K_t(tt+1) = max( min(K_t(tt+1), max(glob.Kgrid)), min(glob.Kgrid));
    
    sig_k_t(tt+1) = std(kdist_t(:,tt), L_t(:,tt));

    % Capital path according to PLM
    
    K_plm_t(tt+1) = exp( (1-A_t(tt))*(param.b_KS(1) + param.b_KS(2).*log(K_t(tt))) ...
        +  A_t(tt)*(param.b_KS(3) + param.b_KS(4).*log(K_t(tt))) );
%     K_plm_t(tt+1) = exp( (1-A_t(tt))*(param.b_KS(1) + param.b_KS(2).*log(K_plm_t(tt))) ...
%         +  A_t(tt)*(param.b_KS(3) + param.b_KS(4).*log(K_plm_t(tt))) );  

%     K_plm0_t(tt+1) = exp( (1-A_t(tt))*(param.b_KS0(1) + param.b_KS0(2).*log(K_plm0_t(tt))) ...
%         +  A_t(tt)*(param.b_KS0(3) + param.b_KS0(4).*log(K_plm0_t(tt))) );    
    
   
end


% close all;
figure;
subplot(2,1,1)
plot(K_t(end-100:end))
hold all
plot(K_plm_t(end-100:end))
% hold all
% plot(K_plm0_t(options.burn+1:end))
legend('ALM','PLM','location','southeast')
% xlim([0,2000])
grid on

subplot(2,1,2)
yyaxis left
plot(A_t(end-100:end), 'color','b')
ylim([-0.1, 1.1])
yyaxis right
plot(K_t(end-100:end), 'color','r')
grid on


% subplot(2,1,2)
% cons = reshape(solved.cons, glob.Nk, glob.Ne, glob.NA, glob.NK);
% plot(glob.kgrid, cons(:,:,1,1))
% hold all
% plot(glob.kgrid, cons(:,:,2,1))
% grid on

% subplot(2,1,2)
% L = reshape(L_t(:,1:end-1), glob.Nkf, glob.Ne, options.T);
% plot(glob.kgridf, L(:,2,1000), '-o')
% hold all
% plot(glob.kgridf, L(:,2,1100), '-o')
% hold all
% plot(glob.kgridf, L(:,2,1200), '-o')
% grid on
% % xlim([0, 100])
% legend('T=1000','T=1100','T=1200','location','northeast')

drawnow


% Pack up output:
sim.A_t = A_t;
sim.K_t = K_t; 
sim.K_plm_t = K_plm_t;
sim.Y_t = Y_t; 
sim.r_t = r_t; 
sim.w_t = w_t; 
sim.L_t = L_t;
sim.kdist_t = kdist_t;
sim.sig_k_t = sig_k_t;








end