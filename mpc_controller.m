function u = mpc_controller(x, xr, X_ref_win, U_ref_win, Qf, Q, R, K_term)

if nargin < 6 || isempty(Q)
    Q = diag([2000 2000 1000000 200 200 200 5 50]);
end
if nargin < 7 || isempty(R)
    R = diag([0.1, 0.1, 0.1, 100, 0.5]);
end
if nargin < 8
    K_term = [];
end

N  = size(U_ref_win, 2);
nu = size(U_ref_win, 1);
nx = size(X_ref_win, 1);

dx0 = x - xr;
dx0(7) = atan2(sin(dx0(7)), cos(dx0(7)));

[F,G] = rolloutPrediction(X_ref_win, U_ref_win);

U_ref_stack = reshape(U_ref_win, [], 1);

[H, f] = build_cost(F, G, dx0, Q, R, Qf, U_ref_stack);

[Aineq, bineq] = build_constraints(dx0, X_ref_win, U_ref_win);

% Polytopic terminal set constraint: u_min <= u_ref_N - K*dx_N <= u_max
% Equivalent linear constraint on DeltaU: [K; -K]*G_N*DU <= rhs
if ~isempty(K_term)
    u_min_t = [-3; -2; -2; -1; -1];
    u_max_t = [ 3;  2;  3;  1;  1];  % u3 upper: 3 m/s^2 gives margin above az_ref~1.62
    u_ref_N = U_ref_win(:, end);
    F_N = F(end-nx+1:end, :);
    G_N = G(end-nx+1:end, :);
    A_tc = [ K_term * G_N; -K_term * G_N];
    b_tc = [(u_ref_N - u_min_t) - K_term * F_N * dx0; ...
            (u_max_t - u_ref_N) + K_term * F_N * dx0];
    n_slack = size(Aineq, 2) - N * nu;
    A_tc = [A_tc, zeros(size(A_tc, 1), n_slack)];
    Aineq = [Aineq; A_tc];
    bineq = [bineq; b_tc];
end

options = optimoptions(@quadprog,'Display','off');
[U,~,exitflag] = quadprog(H, f, Aineq, bineq, [], [], [], [], [], options);

if exitflag <= 0
    warning('MPC QP failed, applying reference input');
    U = [zeros(N*nu, 1); zeros(N, 1)];
end

u = U_ref_win(:,1) + U(1:nu);

end
