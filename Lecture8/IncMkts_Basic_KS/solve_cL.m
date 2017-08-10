function eq  = solve_coeff(Y,param,glob,options)
%SOLVE_CL Solve for value function coefficients and stationary distribution  
%-------------------------------------------------
%   Solves for the value function coeffivient vectors (cK,cC,cE) and the
%   stationary distribution matrix L. Solution is conditional on a
%   conjectured level of equilibrium output, Y. 
%
%   INPUTS
%   - Y         = conjectured value of output, Y
%   - P         = conjectured value of aggregate price level, P
%   - param     = parameters 
%   - glob      = includes state space, function space, approximating functions etc
%   - options   = 
%   OUTPUT
%   - eq        = 
%-------------------------------------------------


%% A. Globals 
s           = glob.s; 
sf          = glob.sf;
pPgrid      = glob.pPgrid;

%% Initialise guesses
cKold       = zeros(glob.Ns,1);
cCold       = zeros(glob.Ns,1);
cEold       = zeros(glob.Ns,1);
cold        = [cKold;cCold;cEold];

% Check if previous solution exists
if exist('glob.c','var')
        cold = glob.c;
end

totaltic    = tic;
%% Bellman iteration
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfunc(cold,s,Y,param,glob,options); 
    % 2. Update c
    cK          = glob.Phi\full(v.vK);      % Note: 'full' re-fills a sparse matrix for computations
    cC          = glob.Phi\full(v.vC);
    cE          = glob.Phi\full(v.vE);    
    c           = [cK;cC;cE];
    % 3. Compute distance and update
    dc          = norm(c-cold)/norm(cold); 
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
end

%% Newton iterations
if strcmp(options.print,'Y');
    fprintf('~~~~~ Newton iterations ~~~~~\n');
end
eq.flag.cconv = false;
for citer = (1:options.Nnewt)
    % 1. Compute values
    [v,jac]     = solve_valfunc_noagg(cold,s,Y,param,glob,options);
    % 2. Update c 
    cKold       = cold(1:glob.Ns); 
    cCold       = cold(glob.Ns+1:2*glob.Ns);
    cEold       = cold(2*glob.Ns+1:end);
    c           = cold - jac\([glob.Phi*cKold - full(v.vK) ;
                               glob.Phi*cCold - full(v.vC) ;
                               glob.Phi*cEold - full(v.vE)]);  
    % 3. Compute distances and update
    dc          = norm(c-cold)/norm(cold);
    cold        = c;
    if strcmp(options.print,'Y');
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
    % 4. Check convergence
    if (dc<options.tolc)
        eq.flag.cconv = true;
    end
    if eq.flag.cconv
        break
    end
end



%% Pack-up output
eq.v    = v;
eq.c    = c;



end

