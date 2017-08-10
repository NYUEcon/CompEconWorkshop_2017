%------------------------------
%   Solve for the policy functions and stationary equilibrium of the 
%   consumption savings problem with incomplete markets, as in 
%   Den Haan, Judd, Juillard (JEDC, 2010). 
%
%   This model attempts to make Krusell Smith work. (Wish me luck...).
%
%   It works!!!1!Shift!1!
%
%   Uses the CompEcon package by Miranda and Fackler.
%
%   James Graham
%   22 July 2016
%
%--------------------------------

close all
clear
clc

%% Set options

% Model options 
options.discmethod     = 'T';      % R=Rouwenhorst, T=two states

% Tolerances, iterations
options.Nbell       = 10;        % Number of Bellman (Contraction) iterations
options.Nnewt       = 30;       % Maximum number of Newton steps
options.tolc        = 1e-8;     % Tolerance on value functions
options.tolgolden   = 1e-8;     % Tolerance for golden search
options.tolL        = 1e-8;    % Tolerance for solving stationary dist
options.itermaxL    = 5000;     % Maximum iterations to find stationary dist L
options.tolKeq      = 1e-4;     % tolerance for solving equilibrium
options.itermaxKeq  = 1000;     % Maximum iterations to find equilibrium

% Tolerances for KS algorithm...
options.T           = 1500;
options.burn        = 500;
options.tolKS       = 1e-4;

% Print / plot 
options.print       = 'Y';      % Print out c-solution convergence
options.eqmprint    = 'Y';      % Print out c-solution convergence
options.fontsize    = 14;       % Plot fontsize
options.linesize    = 1;       
options.plotpolicyfun = 'Y';
options.plotstatdist  = 'Y';

%% Model parameters
param.beta      = 0.99;     % Rate of time preference
param.sigma     = 1;        % Risk aversion parameter
param.alpha     = 0.36;     % Weight on capital in prodution
param.delta     = 0.025;    % depreciation
param.lbar      = 1/0.9;    % units of labor provided
param.Da        = 0.01;     % Size of aggregate fluctuation in productivity
param.mu        = 0.15;     % unemployment benefits
param.u_b       = 0.1;      % unemployment rate in bad state
param.u_g       = 0.04;      % unemployment rate in good state
param.Lbar_b    = 1 - param.u_b;  % employment rate in bad state
param.Lbar_g    = 1 - param.u_g;  % employment rate in good state
% param.Lbar      = 1 - param.ubar;      % aggregate employed

% Steady state values when a is in bad state 
param.Kss       = ((1/param.beta-(1-param.delta))/param.alpha)^(1/(param.alpha-1));
param.rss       = 1/param.beta -1 + param.delta;

% Initial guess for Krusell Smith coefficients 
% param.b_KS0 = [0.123815, 0.965565, 0.137800, 0.963238];
param.b_KS = [0, 1, 0, 1];
param.b_KS = param.b_KS0;

%% Statespace parameters:   [Nk, Ne, NA, NK]
glob.n          = [150, 2, 2, 4];   % Number of nodes in each dimension: 
glob.nf         = [1000, 2, 2, 4];   % Number of points for histogram in stationary distribution
glob.spliorder  = [3, 1, 1, 1];    % Order of splines (use linear if exogenous vars are discrete (not AR1))

glob.kmin       = 0;               % Lower bound on asset holdings
glob.kmax       = 1000;            % Upper bound on asset holdings (suggested by MMV)
glob.kmaxf      = 500;

glob.Kmin       = 35;              % Lower bound on asset holdings
glob.Kmax       = 45;              % Upper bound on asset holdings

%% Setup: basis matrices, expectations, etc

glob = setup(param,glob,options);

%% Solve Krusell Smith algorithm

% Check if previous solution exists
if exist('KS_coeffs.mat','file') == 2
    load('KS_coeffs.mat')
    param.b_KS = b_KS;
end

[ b_KS, R2bad, R2good ] = solve_krusellsmith(param,glob,options);
param.b_KS = b_KS;
save('KS_coeffs.mat','b_KS')   % save coefficients to file



%% Impulse response function
% Warning: need a very large number of burn in periods holding the
% aggregate state at A=0 in order to get the distribution to settle down.
% Hence, the IRF takes a long time.
options.T           = 20000;
options.burn        = 19000;
prevdate=20;
enddate=150;
irf_sim = irf(prevdate, enddate, param, glob, options);

% Save figure to PDF
h = gcf;
set(h,'PaperOrientation','landscape');
print('IRF', '-dpdf', '-bestfit')

%% Get moments from a long simulation

options.T           = 12000;
options.burn        = 2000;

sim = simulation(param,glob,options);

% Compute errors 
errors = (sim.K_t - sim.K_plm_t)./sim.K_plm_t;
fprintf('Max error = %0.4f \n', max(errors*100))
fprintf('Mean error = %0.4f \n', mean(errors*100))

% Compute business cycle statistics
corr.y_sig = corrcoef(sim.Y_t(options.burn+1:options.T), sim.sig_k_t(options.burn+1:options.T));
corr.y_K = corrcoef(sim.Y_t(options.burn+1:options.T),sim.K_t(options.burn+1:options.T));

enddate=150;
figure
subplot(2,2,1)
plot(sim.A_t(options.burn+1:options.burn+enddate))
ylim([-0.1, 1.1])
title('TFP')
grid on

subplot(2,2,2)
plot(sim.K_t(options.burn+1:options.burn+enddate))
grid on
title('Aggregate capital stock')

subplot(2,2,3)
plot(sim.Y_t(options.burn+1:options.burn+enddate))
grid on
title('Aggregate output')

subplot(2,2,4)
plot(sim.sig_k_t(options.burn+1:options.burn+enddate))
title('Std. dev. of wealth')
grid on





%% Solve for a given value of aggregate capital 

% Solve coefficients
solved = solve_coeff(param,glob,options);

options.linesize = 2;

figure
subplot(1,2,1)
kstar = reshape(solved.kstar, length(glob.kgrid), glob.Ne, glob.NA, glob.NK);
plot(glob.kgrid, kstar(:,1,1,2), 'linewidth', options.linesize)
hold all
plot(glob.kgrid, kstar(:,2,1,2), 'linewidth', options.linesize)
hold all
plot(glob.kgrid, kstar(:,1,2,2), 'linewidth', options.linesize)
hold all
plot(glob.kgrid, kstar(:,2,2,2), 'linewidth', options.linesize)
line([0 glob.kmax],[0 glob.kmax], 'color','k', 'linestyle','--')
title('ASSET FUNCTION','fontsize',options.fontsize)
xlabel('Assets today','fontsize',options.fontsize)
ylabel('Capital choice','fontsize',options.fontsize)
legend('(Unemp, Bad)', '(Emp, Bad)','(Unemp, Good)','(Emp, Good)','Location','northwest')
set(gca, 'fontsize', options.fontsize)
grid minor
xlim([0, 40])

subplot(1,2,2)
cons = reshape(solved.cons, length(glob.kgrid), glob.Ne, glob.NA, glob.NK);
plot(glob.kgrid, cons(:,:,1,1), 'linewidth', options.linesize)
hold all
plot(glob.kgrid, cons(:,:,2,1), 'linewidth', options.linesize)
title('CONSUMPTION FUNCTION','fontsize',options.fontsize)
xlabel('Assets today','fontsize',options.fontsize)
ylabel('Consumption choice','fontsize',options.fontsize)
set(gca, 'fontsize', options.fontsize)
grid minor
xlim([0, 40])


% Save figure to PDF
h = gcf;
set(h,'PaperOrientation','landscape');
print('decision_rules', '-dpdf', '-bestfit')


%% Solve for a given value of aggregate capital 

% Solve coefficients
solved = solve_coeff(param,glob,options);

% Simulation
options.T           = 12000;
options.burn        = 2000;
sim = simulation(param,glob,options);

T = 2659;
K = sim.K_t(T);
A = sim.A_t(T);
L = reshape(sim.L_t(:,T), glob.Nkf, glob.Ne);



%% 
options.linesize = 2

figure
subplot(1,2,1)
kstar = reshape(solved.kstar, length(glob.kgrid), glob.Ne, glob.NA, glob.NK);
yyaxis left
plot(glob.kgrid, kstar(:,1,2,2), 'linewidth', options.linesize, ...
                                            'linestyle','-','color','red')
hold all
plot(glob.kgrid, kstar(:,2,2,2), 'linewidth', options.linesize, ...
                                            'linestyle','-','color','blue')
line([0 glob.kmax],[0 glob.kmax], 'color','k', 'linestyle','--')

yyaxis right
plot(glob.kgridf, L(:,2),'-o', 'color','blue')
hold all
plot(glob.kgridf, L(:,1),'-o', 'color','red')
ylim([0, 0.03])

title('ASSET FUNCTION','fontsize',options.fontsize)
xlabel('Assets today','fontsize',options.fontsize)
ylabel('Capital choice','fontsize',options.fontsize)
legend('(Unemp, Good)','(Emp, Good)','Location','northwest')
set(gca, 'fontsize', options.fontsize)
grid minor
xlim([0, 40])

subplot(1,2,2)
cons = reshape(solved.cons, length(glob.kgrid), glob.Ne, glob.NA, glob.NK);
yyaxis left
plot(glob.kgrid, cons(:,1,1,2), 'linewidth', options.linesize, ...
                                            'linestyle','-','color','red')
hold all
plot(glob.kgrid, cons(:,2,2,2), 'linewidth', options.linesize, ...
                                            'linestyle','-','color','blue')

yyaxis right
plot(glob.kgridf, L(:,2),'-o', 'color','blue')
hold all
plot(glob.kgridf, L(:,1),'-o', 'color','red')
ylim([0, 0.03])

title('CONSUMPTION FUNCTION','fontsize',options.fontsize)
xlabel('Assets today','fontsize',options.fontsize)
ylabel('Consumption choice','fontsize',options.fontsize)
legend('(Unemp, Good)','(Emp, Good)','Location','northwest')
set(gca, 'fontsize', options.fontsize)
grid minor
xlim([0, 40])


% Save figure to PDF
h = gcf;
set(h,'PaperOrientation','landscape');
print('decision_rule_hist', '-dpdf', '-bestfit')




