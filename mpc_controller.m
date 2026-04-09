function u = mpc_controller(x, xr, X_ref_win, U_ref_win, Qf, Q, R)

if nargin < 6 || isempty(Q)
    Q = diag([2000 2000 1000000 200 200 200 5 50]);
end
if nargin < 7 || isempty(R)
    R = diag([0.1, 0.1, 0.1, 100, 0.5]);
end

N  = size(U_ref_win, 2);
nu = size(U_ref_win, 1);

dx0 = x - xr;
% Wrapping
dx0(7) = atan2(sin(dx0(7)), cos(dx0(7)));

[F,G] = rolloutPrediction(X_ref_win, U_ref_win);

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
    U = [zeros(N*nu, 1); zeros(N, 1)];
end

u = U_ref_win(:,1) + U(1:nu);

end
