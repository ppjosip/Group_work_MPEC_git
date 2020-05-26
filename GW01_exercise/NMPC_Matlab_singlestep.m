function  out = NMPC_Matlab_singlestep(qp_W,qp_gradJ,qp_gradhT,qp_h,FJk,parNMPC,parSim);

% Initial guess

% Initialization
optVars   	= parNMPC.optVars0;
lambda      = parNMPC.lambda0;      % Dual variables for linear inequality
                                    % constraints (as used in qpOASES):
                                    % lbA <= Ax <= ubA
mu          = parNMPC.mu0;          % Dual variables for lower and upper 
                                    % bounds on optimization variables:
                                    % lb <= x <= ub
kktCond     = Inf*ones(parNMPC.nOptVars+parNMPC.nConstr,1);
kSQP     	= 1;

% Pre-define
p = [parNMPC.x0; parNMPC.uprev; parSim.ref.p_im.data(1); parSim.ref.x_bg.data(1); parNMPC.Q1_val; parNMPC.Q2_val; parNMPC.R1_val; parNMPC.R2_val];

% SQP iterations until stopping criteria reached
tic
while (norm(kktCond,2) >= parNMPC.kktTol) && (kSQP <= parNMPC.nSQP)
    
    % Get numerical values for QP matrices/vectors and bounds
    W       = full(qp_W(optVars,[p; lambda]));
    gradJ 	= full(qp_gradJ(optVars, p));
    gradhT  = full(qp_gradhT(optVars, p));
    h     	= full(qp_h(optVars, p));
    
    % Set proper bounds for s_xk = x(k+1) - x(k)
    lb_sxk 	= parNMPC.lbx - optVars;
    ub_sxk	= parNMPC.ubx - optVars;
    
    % Evaluation of KKT conditions
    kktCond = [gradJ - gradhT'*lambda - mu; h];
    
    % Solve the QP
    if kSQP==1
        [QP,sQP,~,exitflag,~,lambdaQP,~] = qpOASES_sequence('i',W,gradJ,...
            gradhT,lb_sxk,ub_sxk,-h,-h);
    else
        [sQP,~,exitflag,~,lambdaQP,~] = qpOASES_sequence('m',QP,W,gradJ,...
            gradhT,lb_sxk,ub_sxk,-h,-h);
    end % if
    lambda_OptVars  = lambdaQP(1:length(lb_sxk));
    lambda_AOptVars = lambdaQP(length(lb_sxk)+1:end);
    
    % Update optimisation variables and lambda/mu
    xprev   = optVars;          % for debugging
    sx      = sQP;
    smu     = lambda_OptVars - mu;
    slambda = lambda_AOptVars - lambda;
    optVars = optVars + parNMPC.alphaSQP*sx;
    mu      = mu + parNMPC.alphaSQP*smu;
    lambda  = lambda + parNMPC.alphaSQP*slambda;
    
    % Figure during loop (for debugging)
    if 0
        figure(99); clf;
        ax(1) = subplot(2,1,1); hold on; box on; grid on;
        plot([0 nOptVars],-0.25*[1 1],'b--');
        plot([0:3:nOptVars],[NaN; lb(2:3:end)],'b-');
        plot([0:3:nOptVars],[options.x0(1); xprev(2:3:end)],'bx');
        plot([0:3:nOptVars],[options.x0(1); optVars(2:3:end)],'b.','MarkerSize',10);
        plot([0:3:nOptVars],[options.x0(2); xprev(3:3:end)],'rx');
        plot([0:3:nOptVars],[options.x0(2); optVars(3:3:end)],'r.','MarkerSize',10);
        plot([0:3:nOptVars],[NaN; sQP(3:3:end)],'ko','LineWidth',1);
        ylabel('$x$');
        ax(2) = subplot(2,1,2); hold on; box on; grid on;
        plot([0 nOptVars],-[1 1],'k--');
        plot([0 nOptVars],[1 1],'k--');
        plot([0:3:nOptVars],[lb(1:3:end); NaN],'k-');
        plot([0:3:nOptVars],[ub(1:3:end); NaN],'k-');
        plot([0:3:nOptVars],[xprev(1:3:end); NaN],'kx');
        plot([0:3:nOptVars],[optVars(1:3:end); NaN],'k.','MarkerSize',10);
        ylabel('$u$');
        xlabel('Index $k$');
        disp(exitflag);
        pause(0.5);
    end % if
    
    % Update iteration index
    kSQP = kSQP + 1;
    
end % while
time_to_solve = toc;

%Generate Solution
out.optVarsPred = optVars;

end