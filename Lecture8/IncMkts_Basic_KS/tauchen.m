function [T, grid] = tauchen(N,rho,sigma_eps,x_bar)
%TAUCHEN - create transition matrix for an AR(1) using the 
%          Tauchen discretization method
%----------------------------------------------------------------
%   The continuous AR(1) process is given by: 
%   x_t = x_bar + rho*x_t-1 + eps_t 
%   This process is then descretized using Tauchen method.
%
%   [T, grid] = tauchen(N,rho,sigma_eps,x_bar)
%
%   INPUTS:
%       N = number of grid points/states
%       rho = AR(1) coefficient
%       sigma_eps = standard deviation of error term
%       x_bar = constant term
%   
%   OUTPUTS:
%       T = transition matrix
%       grid = a vector containing the state space
%------------------------------------------------------------------
    
% Set Tauchen model parameters
x_var =  sigma_eps^2/(1-rho^2);
x_mean = x_bar/(1-rho);
std_x = sqrt(x_var);                      % create std of the x process
x_N = x_mean + 3*std_x;
x_1 = x_mean - 3*std_x;
grid = linspace(x_1, x_N, N);            % equispaced grid between the two boundary points
d = grid(2)-grid(1);                            % distance between two successive points
  
% Create transition matrix
T = zeros(N,N);
for j = 1:N     % moving from state j, to state k
    T(j,N) = 1 - normcdf((grid(N) - d/2 - x_bar - rho*grid(j))/sigma_eps, 0, 1);
    T(j,1) = normcdf((grid(1) + d/2 - x_bar - rho*grid(j))/sigma_eps, 0, 1);
    for k = 2:N-1
        T(j,k) = normcdf((grid(k) + d/2 - x_bar - rho*grid(j))/sigma_eps, 0, 1)...
            - normcdf((grid(k) - d/2 - x_bar - rho*grid(j))/sigma_eps, 0 , 1);
    end
end

% Adjusts row probabilities to ensure that they sum to one 
for i = 1:N
    T(i,:) = T(i,:)/sum(T(i,:));   
end





end

