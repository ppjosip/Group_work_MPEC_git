function [qp_W,qp_gradJ,qp_gradhT,qp_h,FJk] = createCasadiFunctions(parNMPC);
import casadi.*

%% Defintion of model dynamics and objective funtions
% States and inputs
x       = MX.sym('x', parNMPC.nStates, 1);
u       = MX.sym('u', parNMPC.nInputs, 1);

% System dynamics
[xdot, y] = ROModel(x, u);

% Objective function 
r       = MX.sym('r', 2, 1);
du      = MX.sym('du', parNMPC.nInputs, 1);
Q1      = MX.sym('Q1');
Q2      = MX.sym('Q2');
R1      = MX.sym('R1');
R2      = MX.sym('R2');
J      = Q1*(y(1)-r(1))^2 + Q2*(y(2)-r(2))^2 + R1*du(1)^2 + R2*du(2)^2;

% Create CasADI function
fxdot = Function ('fxdot', {x,u}, {xdot});
fJDisc    = Function ('fJDisc', {x,du,r,[Q1;Q2;R1;R2]},{J});

%% Integration/discretization using RK4
xStart = MX.sym('xStart', parNMPC.nStates, 1);
u      = MX.sym('u', parNMPC.nInputs, 1);
TRK4   = parNMPC.Ts/parNMPC.nRK4;

% Loop over intervals for fourth order Runge-Kutta discretization methode
xEnd = xStart;
JEnd = 0;

for l = 1:parNMPC.nRK4
    k1x = fxdot(xEnd, u);
    k2x = fxdot(xEnd + TRK4/2 * k1x, u);
    k3x = fxdot(xEnd + TRK4/2 * k2x, u);
    k4x = fxdot(xEnd + TRK4 * k3x, u);
    xEnd = xEnd + TRK4/6 * (k1x + 2*k2x + 2*k3x + k4x);
end

% Create CasADi functions
fxDisc  = Function('fxDisc', {xStart,u}, {xEnd});

%% Construct NLP
% Intialization
optVars     = [];   % Vector of optimization variables (i.e. states and inputs)
optVars0    = [];   % Initial guess for optimization variables
lb          = [];   % Lower bound of optimization variables
ub          = [];   % Upper bound of optimization variables
Jk          =  0;   % Initialization of objective function
h           = [];   % (In-)equality constraints

% Pre-define CasADi variables
U = MX.sym('U_', parNMPC.nInputs, 1, parNMPC.N);
S = MX.sym('S_', parNMPC.nStates, 1, parNMPC.N+1);
% Pre-define parameters - p = [x0; uprev; r; Q1; Q2; R1; R2]
p = MX.sym('p', parNMPC.nStates + parNMPC.nInputs + parNMPC.nOutputs+4); 

% Construct NLP step-by-step
for k = 1:parNMPC.N+1
    
    % System dynamics and objective function
    if k==parNMPC.N+1
        % Skip, because at final time step, no further integration of
        % system dynamics necessary.
    else
        % Integration of system dynamics and objective function
        if k==1
            % Hardcoded initial condition
            XEnd = fxDisc(p(1:5),U{k});
            Jk = fJDisc(p(1:5), U{k}-p(6:7), p(8:9), p(10:13));
        else
            XEnd = fxDisc(S{k},U{k});
            Jk = Jk + fJDisc(S{k},U{k}-U{k-1}, p(8:9), p(10:13));
        end % if
        
        % Add equality constraint for continuity (i.e. closing gaps):
        h = [h;  XEnd - S{k+1}];
        
    end % if
    
    % States
    if k==1
        % Skip, because we have hardcoded the initial condition above
    else
        % Add states to vector of optimization variables
        optVars = [optVars; S{k}];
        
        % Lower- and upper bounds for states
        lb = [lb; -inf; -inf; -inf; -inf; -inf];
        ub = [ub; inf; inf; inf; inf; inf]; % No State constraints
        
        % Add initial guess of states
        optVars0 = [optVars0; p(1:5)];
        
    end % if
    
    % Inputs
    if k==parNMPC.N+1
        % Skip, because no control input at final time step
    else
        % Add inputs to vector of optimization variables
        optVars = [optVars; U{k}];
        
        % Lower- and upper bounds for inputs
        lb = [lb; 0; 0];
        ub = [ub; 1; 1]; 
        
        % Add initial guess for inputs
        optVars0 = [optVars0; 0; 0];
        
    end % if
    
end % for

%% SQP
% Construct Lagrange function
lambda = MX.sym('lambda', numel(h));
L = Jk + lambda'*h;

W = hessian(L, optVars);
gradJ = jacobian(Jk, optVars)';
gradhT = jacobian(h, optVars);

qp_W = Function('qp_W', {optVars, [p;lambda]}, {W});
qp_gradJ = Function ('qp_gradJ', {optVars, p}, {gradJ});
qp_gradhT = Function ('qp_gradhT', {optVars, p}, {gradhT});
qp_h = Function ('qp_h', {optVars, p}, {h});

%FJK
FJk = 0;

end