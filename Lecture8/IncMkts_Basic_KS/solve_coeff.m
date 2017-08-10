function solved  = solve_coeff(param,glob,options)
%SOLVE_CL Solve for value function coefficients and stationary distribution  
%-------------------------------------------------
%   Solves for the value function coeffivient vectors (cC,cE). 
%
%   INPUTS
%   - param     = parameters 
%   - glob      = includes state space, function space, approximating functions etc
%   - options   = 
%   OUTPUT
%   - solved        = 
%-------------------------------------------------


%% A. Globals 
s          = glob.s; 

%% Initialise guesses
cVold       = zeros(glob.Ns,1);
cEold       = zeros(glob.Ns,1);
cold        = [cVold;cEold];

% Check if previous solution exists
if exist('vfi_coeffs.mat','file') == 2
    load('vfi_coeffs.mat')
    % Check if coefficient vector is the correct length
    if (size(c,1) == 2*glob.Ns)
        cold = c;
    end
end


% Check if previous solution exists
% if exist('glob.c','var')
%         cold = glob.c;
% end

totaltic    = tic;
%% Bellman iteration
if strcmp(options.print,'Y')
    fprintf('~~~~~ Bellman iterations ~~~~~\n');
end
for citer = (1:options.Nbell)
    glob.citer  = citer;
    % 1. Compute values;
    v           = solve_valfunc(cold,s,param,glob); 
    % 2. Update c
    cV          = glob.Phi\full(v.vV);      % Note: 'full' re-fills a sparse matrix for computations
    cE         = glob.Phi\full(v.vE);
    c           = [cV; cE];
    % 3. Compute distance and update
    dc          = max(abs(c-cold));     
%     dc          = norm(c-cold)/norm(cold); 
    cold        = c;
    if strcmp(options.print,'Y')
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
end

%% Newton iterations
if strcmp(options.print,'Y')
    fprintf('~~~~~ Newton iterations ~~~~~\n');
end
solved.flag.cconv = false;
for citer = (1:options.Nnewt)
    % 1. Compute values
    [v,jac]     = solve_valfunc(cold,s,param,glob);
    % 2. Update c 
    cVold       = cold(1:end/2); 
    cEold      = cold(end/2+1:end);
    c           = cold - jac\([glob.Phi*cVold - full(v.vV) ;
                               glob.Phi*cEold - full(v.vE) ]);  
    % 3. Compute distances and update
    dc          = max(abs(c-cold));
    cold        = c;
    if strcmp(options.print,'Y')
        fprintf('%i\tdc = %1.2e\tTime: %3.2f\n',citer,dc,toc(totaltic));
    end
    % 4. Check convergence
    if (dc<options.tolc)
        solved.flag.cconv = true;
    end
    if solved.flag.cconv
        fprintf('\n Solved individual problem \n');
        fprintf('-----------------\n')
        break
    end
end

%% Find coefficients for kprime
ck          = glob.Phi\v.kstar;

%% Create consumption function
cons = menufun('consumption',s,v.kstar,param,glob);
cc   = glob.Phi\cons;


%% Pack-up output
solved.v    = v;
solved.c    = c;
solved.ck   = ck;
solved.cc   = cc;
solved.kstar = v.kstar;
solved.cons = cons;


%% Save value function coefficients to file
solved.coeffs.c    = c;
save('vfi_coeffs.mat','c')


end

