function [glob] = setup(param,glob,options)
%SETUP Prepares model objects for collocation procedure, no agg uncertainty
%-------------------------------------------------
%   This file prepares the Markov process (continuous or discretized), the
%   state space (across exogenous and endogenous variables), the function
%   space for approximating with, elements for use in approximation of the
%   stationary distribution, and basis functions for use in the collocation
%   procedure computations
%
%   INPUTS
%   - param     = model parameters
%   - glob      = global variables, including those used in the approximating functions
%   - options   = options, including continuous vs. discrete Markov, Tauchen vs. Rouwenhurst, 
%-------------------------------------------------
fprintf('-----------------\n');
fprintf('Setup\n');

%% State space for idiosyncratic productivity
% One persistent shock
Nk = glob.n(1);
Ne = glob.n(2);
NA = glob.n(3);
NK = glob.n(4);

switch options.discmethod
    case 'R'
        [P,egrid,Psse]      = setup_Markov(Ne,param.sig_e,param.rho_e,1);
    case 'T'
        glob.P = [0.525, 0.35, 0.03125, 0.09375; ...  
             0.038889, 0.836111, 0.002083, 0.122917; ...
             0.09375, 0.03125, 0.291667, 0.583333; ...
             0.009115, 0.115885, 0.024306, 0.850694];
        Psse = glob.P^1000;
        Psse = Psse(1,:);
        
        egrid = [0, 1];        
        Agrid = [0, 1];
        
%         P               = [0.6, 0.4; 0.044445, 0.955555];
%         Psse = P^1000;
%         Psse = Psse(1,:);
%         egrid           = [0, 1];
end

egrid0          = egrid;
Agrid0          = Agrid;

%% Derived transition matrices

% Transition matrix for aggregate state only
glob.PA = zeros(2,2);  
glob.PA(1,1) = glob.P(1,1)+glob.P(1,2); 
glob.PA(1,2) = 1 - glob.PA(1,1);  
glob.PA(2,2) = glob.P(3,3) + glob.P(3,4); 
glob.PA(2,1) = 1 - glob.PA(2,2);

% Conditional employment transition matrices
glob.Pe_b = zeros(2,2);
glob.Pe_b(1,1) = glob.P(1,1)+glob.P(1,3);
glob.Pe_b(1,2) = 1 - glob.Pe_b(1,1);
glob.Pe_b(2,2) = glob.P(2,2) + glob.P(2,4); 
glob.Pe_b(2,1) = 1 - glob.Pe_b(2,2);

glob.Pe_g = zeros(2,2);
glob.Pe_g(1,1) = glob.P(3,1)+glob.P(3,3);
glob.Pe_g(1,2) = 1 - glob.Pe_g(1,1);
glob.Pe_g(2,2) = glob.P(4,2) + glob.P(4,4); 
glob.Pe_g(2,1) = 1 - glob.Pe_g(2,2);

%% Set up grid-space for endogenous variables
spliorder       = glob.spliorder;
scale          = 0.5; 

fprintf('Using a polynomially spaced grid for `k`\n');
x=linspace(0,0.5,Nk)';
y=x.^7/max(x.^7);
kgrid=glob.kmin + (glob.kmax - glob.kmin)*y;

% fprintf('Using a logged grid space\n');
% kgrid           = exp(nodeunif(Nk,log(glob.kmin+scale),...
%                         log(glob.kmax+scale) ))-scale;  % Adds curvature

kgrid0          = kgrid;    % Save for computing basis matrices in valfunc.m:line9   

%% Set up grid for aggregate capital 
fprintf('Using an evenly spaced grid for `K`\n');
Kgrid           = nodeunif(NK, glob.Kmin, glob.Kmax);  % Adds curvature                    
Kgrid0          = Kgrid;

%% Function space and nodes (fspace adds knot points for cubic splines)
fspace          = fundef({'spli',kgrid,0,spliorder(1)},...
                         {'spli',egrid,0,spliorder(2)},...
                         {'spli',Agrid,0,spliorder(3)},...
                         {'spli',Kgrid,0,spliorder(4)});
sgrid           = funnode(fspace);
s               = gridmake(sgrid);
Ns              = size(s,1);

%% Reconstruct grids after fspace added points for the spline (adds two knot points for cubic spline)
kgrid           = sgrid{1};  
egrid           = sgrid{2};  
Agrid           = sgrid{3};  
Kgrid           = sgrid{4}; 

Nk              = size(kgrid,1);
Ne              = size(egrid,1); 
NA              = size(Agrid,1);
NK              = size(Kgrid,1);

%% Set up fine grid-space for endogenous variables 
Nkf             = glob.nf(1);

% fprintf('Using a logged grid space\n');
% kgridf           = exp(nodeunif(Nkf,log(glob.kmin+scale),...
%                         log(glob.kmax+scale) ))-scale;

fprintf('Using an evenly spaced fine grid for `k`\n');
kgridf          =   nodeunif(Nkf,glob.kmin,glob.kmaxf);                    
egridf          = egrid;
Agridf          = Agrid;
Kgridf          = Kgrid;

sf              = gridmake(kgridf,egridf,Agridf,Kgridf);
Nsf             = size(sf,1);

Nkf             = size(kgridf,1);
Nef             = size(egridf,1);
NAf             = size(Agridf,1);
NKf             = size(Kgridf,1);


%% Pre-Compute expectations matrix
glob.Phi        = funbas(fspace,s);
glob.Phif       = funbas(fspace, sf);

%% Compute QA matrix for approximation of stationary distribution
% FIX THIS>>>

glob.QE         = kron(glob.P,ones(Nkf,1)); 

%% Create basis matrices for expectations matrix (need only compute once)
Phi_k      = splibas(kgrid0,0,spliorder(1),s(:,1));
Phi_e      = splibas(egrid0,0,spliorder(2),s(:,2));
Phi_A      = splibas(Agrid0,0,spliorder(3),s(:,3));
glob.Phi_K      = splibas(Kgrid0,0,spliorder(4),s(:,4));

glob.Phi_eA    = dprod(Phi_A, Phi_e);
glob.Phi_eAK   = dprod(glob.Phi_K, glob.Phi_eA);
glob.Phi_keA   = dprod(Phi_A, dprod(Phi_e, Phi_k));

glob.Emat       = kron(speye(NK), kron(glob.P,speye(Nk)));
% spy(glob.Emat)    % Can let us look at these sparse matrices


% Compute for the fine grid 
% Phi_ef          = splibas(egrid0,0,spliorder(2),sf(:,2));
% Phi_Af          = splibas(Agrid0,0,spliorder(3),sf(:,3));
% Kprimef         = menufun('Kprime', sf, [], param, glob);
% Phi_Kf          = splibas(Kgrid0,0,spliorder(4),Kprimef);
% 
% glob.Phi_eAKf = dprod(Phi_Kf, dprod(Phi_Af, Phi_ef));


% glob.Phi_ef     = splibas(egrid0,0,spliorder(2),sf(:,2));       % Used when solving on fine grid

%% Declare additional global variables
glob.kgrid0    = kgrid0;          % unique elements of pP grid
glob.kgrid     = kgrid;           % full state space for pP grid
glob.kgridf     = kgridf;
glob.egrid0     = egrid0;           % unique elements of a grid
glob.egrid      = egrid;            % full state space for a grid
glob.egridf     = egridf;
glob.Agrid      = Agrid;
glob.Agrid0     = Agrid0;
glob.Agridf     = Agridf;
glob.Kgrid      = Kgrid;
glob.Kgrid0     = Kgrid0;
glob.Kgridf     = Kgridf;
glob.Psse       = Psse;             % stationary distribution for employment
glob.Ne         = Ne;               % length of state space grid for a 
glob.Nk         = Nk;              
glob.NA         = NA;              
glob.NK         = NK;              
glob.Nkf        = Nkf;              % length of state space grid for pP
glob.fspace     = fspace;           % function space object for the model
glob.s          = s;                % full state space 
glob.sf         = sf;
glob.Ns         = Ns;               % size of full state space
glob.Nsf        = Nsf;

fprintf('Setup complete\n');
fprintf('-----------------\n');


end


        
   
