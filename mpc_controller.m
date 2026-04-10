function u = mpc_controller(x, xr, X_ref_win, U_ref_win, Qf)

% N  = size(U_ref_win,2);
nu = size(U_ref_win,1);

dx0 = x - xr;

% Rollout prediction
[F,G] = rolloutPrediction(X_ref_win, U_ref_win);

% Weight matrices
%Q  = diag([10, 10, 1e-6, 1, 1, 0.1, 5, 1]);
%Q = diag([2000 2000 1000000 200 200 200 5 50]);
Q = diag([20000 20000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 100, 0.5]);

% Stack reference input over horizon
U_ref_stack = reshape(U_ref_win, [], 1);

% Build cost
[H, f] = build_cost(F, G, dx0, Q, R, Qf, U_ref_stack);

% Build constraints
[Aineq, bineq] = build_constraints(dx0, X_ref_win, U_ref_win);

% Solve QP
options = optimoptions(@quadprog,'Display','off');
[U,~,exitflag] = quadprog(H, f, Aineq, bineq, [], [], [], [], [], options);

% Fallback if solver fails
if exitflag <= 0
    warning('MPC QP failed, applying reference input');
    U = U_ref_stack;
end

% Apply first control input
u = U(1:nu);

end
