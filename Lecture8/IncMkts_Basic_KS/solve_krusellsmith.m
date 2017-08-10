function [ b_KS, R2bad, R2good ] = solve_krusellsmith(param,glob,options)
%SOLVE_KRUSELLSMITH solves for equilibrium using the KS algorithm
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
damp = 0.75;

b_KS_new = zeros(1,4);

dif_B = 100;
while dif_B >= options.tolKS

% Solve individual problem for current KS parameters param.b_KS
% solved = solve_coeff(param,glob,options);
totaltic    = tic;

sim = simulation(param,glob,options);
fprintf('Simulation time: %3.2f\n',toc(totaltic));
fprintf('-----------------\n')

ibad=0;           % count how many times the aggregate shock was bad
igood=0;          % count how many times the aggregate shock was good
xbad=0;  ybad=0;  % regression-variables for a bad state
xgood=0; ygood=0; % regression-variables for a good state
for tt = options.burn+1:options.T-1
   if sim.A_t(tt)==0
      ibad = ibad + 1;
      xbad(ibad,1) = log(sim.K_t(tt));
      ybad(ibad,1) = log(sim.K_t(tt+1));
   else
      igood = igood + 1;
      xgood(igood,1) = log(sim.K_t(tt));
      ygood(igood,1) = log(sim.K_t(tt+1));
   end
end

% run the OLS regression ln(km')=B(1)+B(2)*ln(km) for a bad agg. state 
% and compute R^2 (which is the first statistic in s5)
[b_KS_new(1:2),~,~,~,s5] = regress(ybad,[ones(ibad,1) xbad]);
R2bad = s5(1);

% run the OLS regression ln(km')=B(3)+B(4)*ln(km) for a good agg. state 
% and compute R^2 (which is the first statistic in s5)
[b_KS_new(3:4),~,~,~,s5] = regress(ygood,[ones(igood,1) xgood]);
R2good = s5(1);

% Check if KS parameters are converging
dif_B = max(abs(param.b_KS - b_KS_new));    % norm(param.b_KS - b_KS_new); 

fprintf('-----------------\n')
fprintf('Norm(KS) = %.8f \n', dif_B)
fprintf('Old KS coeffs: %.4f, %.4f, %.4f, %.4f  \n', ...
    param.b_KS(1), param.b_KS(2), param.b_KS(3), param.b_KS(4))
fprintf('New KS coeffs: %.4f, %.4f, %.4f, %.4f  \n', ...
    b_KS_new(1), b_KS_new(2), b_KS_new(3), b_KS_new(4))
fprintf('-----------------\n')


% Updating via dampening
param.b_KS = damp*param.b_KS + (1-damp)*b_KS_new;
% Update dampening parameter
damp = 0.95*damp;

      

end

% Pack up output 
b_KS = b_KS_new;




end