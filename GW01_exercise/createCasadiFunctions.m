function [qp_W,qp_gradJ,qp_gradhT,qp_h,FJk] = createCasadiFunctions(parNMPC);
import casadi.*

%% Defintion of model dynamics and objective funtions
% States and inputs
x       = MX.sym('x', 2, 1); % 2x1 matrix, state 1 is p_im and state two is x_bg
u       = MX.sym('u', 2, 1); % 2x1 matrix, two inputs --> u_vtg & u_egr
nStates = length(x);
nInputs = length(u);

% System dynamics

% Objective function --> please correct it! 
J       = parNMPC.Q1_val*(x(1)- p_ref)^2 + parNMPC.Q2_val*(x(2) - x_ref)^2 + parNMPC.R1_val*(u(1))^2 + parNMPC.R2_val*(u(2));

% Create CasADI function
fxdot = Function ('fxdot', {x,u}, {xdot});
fJ    = Function ('fJ', {x,u}, {J});

%% Integration/discretization using RK4
xStart = MX.sym('xStart', nStates, 1);
u      = MX.sym('u', nInputs, 1);
TRK4   = parNMPC.Ts/parNMPC.nRK4;

% Loop over intervals for fourth order Runge-Kutta discretization methode
xEnd = xStart;
JEnd = 0;

for l = 1:options.nRK4
    k1x = fxdot(xEnd, u);
    k2x = fxdot(xEnd + TRK4/2 * k1x, u);
    k3x = fxdot(xEnd + TRK4/2 * k2x, u);
    k4x = fxdot(xEnd + TRK4 * k3x, u);
    k1J = fJ(xEnd, u);
    k2J = fJ(xEnd + TRK4/2 * k1x, u);
    k3J = fJ(xEnd + TRK4/2 * k2x, u);
    k4J = fJ(xEnd + TRK4 * k3x, u);
    xEnd = xEnd + TRK4/6 * (k1x + 2*k2x + 2*k3x + k4x);
    JEnd = JEnd + TRK4/6 * (k1J + 2*k2J + 2*k3J + k4J);
end

% Create CasADi functions
fxDisc  = Function('fxDisc', {xStart,u}, {xEnd});
fJDisc  = Function('fJDisc', {xStart,u}, {JEnd});

%% Construct NLP
% Intialization
optVars     = [];   % Vector of optimization variables (i.e. states and inputs)
optVars0    = [];   % Initial guess for optimization variables
lb          = [];   % Lower bound of optimization variables
ub          = [];   % Upper bound of optimization variables
Jk          =  0;   % Initialization of objective function
h           = [];   % (In-)equality constraints

% Pre-define CasADi variables
deltaU = MX.sym('deltaU_', nInputs, 1, parNMPC.N);
X = MX.sym('X_', nStates, 1, parNMPC.N+1);

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
            XEnd = fxDisc(parNMPC.x0, deltaU{k});
            Jk = Jk + fJDisc(parNMPC.x0,deltaU{k});
        else
            XEnd = fxDisc(X{k},deltaU{k});
            Jk = Jk + fJDisc(X{k},deltaU{k});
        end % if
        
        % NEEDS TO BE CORRECTED! compute proper h!
        
        
        % Add equality constraint for continuity (i.e. closing gaps):
        h = [h;  XEnd - X{k+1}];
        
    end % if
    
    % States
    if k==1
        % Skip, because we have hardcoded the initial condition above
    else
        % Add states to vector of optimization variables
        optVars = [optVars; X{k}];
        
        % Lower- and upper bounds for states
        lb = [lb];
        ub = [ub]; % No State constraints
        
        % Add initial guess of states
        optVars0 = [optVars0; options.x0];
        
    end % if
    
    % Inputs
    if k==options.N+1
        % Skip, because no control input at final time step
    else
        % Add inputs to vector of optimization variables
        optVars = [optVars; deltaU{k}];
        
        % Lower- and upper bounds for inputs
        lb = [lb; 0; 0];
        ub = [ub; 1; 1]; 
        
        % Add initial guess for inputs
        optVars0 = [optVars0; 0];
        
    end % if
    
end % for

%% SQP
% Construct Lagrange function
lambda = MX.sym('lambda', numel(h));
L = Jk + lambda'*h;

W = hessian(L, optVars);
gradJ = jacobian(Jk, optVars)';
gradhT = jacobian(h, optVars);

qp_W = Function('qp_W', {optVars, lambda}, {W});
qp_gradJ = Function ('qp_gradJ', {optVars}, {gradJ});
qp_gradhT = Function ('qp_gradhT', {optVars}, {gradhT});
qp_h = Function ('qp_h', {optVars}, {h});

end