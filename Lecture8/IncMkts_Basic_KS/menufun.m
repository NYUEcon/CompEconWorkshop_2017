function out = menufun(flag,s,kprime,param,glob)
%MENUFUN 
%-------------------------------------------------
%
%   INPUTS
%   - flag      = 
%   - s         = state space
%   - x         = 
%   - Y         = conjectured value of output, Y
%   - param     = 
%   - glob      =
%   - options   = 
%   OUTPUT
%   - v         = 
%   - jac       = Jacobian of the value functions  
%-------------------------------------------------

%%
% Globals
    % None
% Parameters
sigma = param.sigma;
delta = param.delta;
alpha = param.alpha;
lbar = param.lbar;
u_b = param.u_b;
u_g = param.u_g;
Lbar_b = param.Lbar_b;
Lbar_g = param.Lbar_g;
mu = param.mu;    
Da = param.Da;
b_KS = param.b_KS;


%% Cases
switch flag
    case 'bounds'
        k = s(:,1);
        e = s(:,2);
        A = s(:,3);
        K = s(:,4);
        r       = interest(A, K);
        w       = wage(A, K);
        tau     = tax(A);
        kmin    = ones(size(s,1),1)*glob.kmin; 
        kmax    = ones(size(s,1),1)*glob.kmax;    % ((1-tau).*lbar.*e + mu.*(1-e)).*w + (1+r-delta).*k; %   % 
        out     = [kmin,kmax];
    case 'felicity'
        k = s(:,1);
        e = s(:,2);
        A = s(:,3);
        K = s(:,4);
        tau     = tax(A);
        r       = interest(A, K);
        w       = wage(A, K);
        c       = ((1-tau).*lbar.*e + mu.*(1-e)).*w + (1+r-delta).*k - kprime;
        out     = utility(c);
    case 'consumption'
        k = s(:,1);
        e = s(:,2);
        A = s(:,3);
        K = s(:,4);
        tau     = tax(A);
        r       = interest(A, K);
        w       = wage(A, K);
        out     = ((1-tau).*lbar.*e + mu.*(1-e)).*w + (1+r-delta).*k - kprime;
    case 'interest'
        A = s(:,3);
        K = s(:,4);
        out      = interest(A, K);
    case 'wage'
        A = s(:,3);
        K = s(:,4);        
        out     = wage(A, K);        
    case 'output'
        A = s(:,3);
        K = s(:,4);        
        out     = output(A, K);    
    case 'Kprime'
        A = s(:,3);
        K = s(:,4);
        out = exp( (1-A).*(b_KS(1) + b_KS(2).*log(K)) ...
                   +  A.*(b_KS(3) + b_KS(4).*log(K)) );
end


%% NESTED FUNCTIONS

% Utility function
function u = utility(c)
    u = zeros(length(c),1);
    if sigma == 1
        % % Log utility, deal with negative consumption
        for ii = 1:length(c)
            if c(ii) <=0
                u(ii,1) = -1e+3 + c(ii);
            else
                u(ii,1) = log(c(ii));
            end
        end
    else
        for ii = 1:length(c)
            if c(ii) <=0
                u(ii,1) = -1e+3 + c(ii);
            else
                u(ii,1) = c(ii).^(1-sigma)/(1-sigma);
            end
        end
    end
end

% Interest rate as function of aggregate capital
function r = interest(A,K)
    prod = (1-A)*(1-Da) + A*(1+Da);
    L = lbar*((1-A)*Lbar_b + A*Lbar_g);
    r = alpha.*prod.*K.^(alpha-1).*(L).^(1-alpha); 
end

% Wage as function of aggregate capital
function w = wage(A,K)
    prod = (1-A)*(1-Da) + A*(1+Da);
    L = lbar*((1-A)*Lbar_b + A*Lbar_g);
    w = (1-alpha).*prod.*K.^(alpha).*(L).^(alpha); 
end

% Output as function of aggregate capital
function y = output(A,K)
    prod = (1-A)*(1-Da) + A*(1+Da);
    L = lbar*((1-A)*Lbar_b + A*Lbar_g);
    y = prod.*K.^alpha.*(L).^(1-alpha);  
end

% taxation as function of aggregate capital
function tau = tax(A)
    unemployed = ((1-A)*(u_b) + A*(u_g));
    L = lbar*((1-A)*Lbar_b + A*Lbar_g);
    tau = mu.*unemployed./(L);
end




end
        
        
        
