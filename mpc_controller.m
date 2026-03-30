function u = mpc_controller(x, xr, X_ref_win, U_ref_win, Qf)

N = 15;

dx0 = x - xr;
% Wrapping yaw error
dx0(7) = atan2(sin(dx0(7)), cos(dx0(7)));

[F,G] = rolloutPrediction(X_ref_win, U_ref_win);

nx = size(Ad,1);
nu = size(Bd,2);

Q  = diag([10, 10, 0.1, 1, 1, 0.1, 5, 1]);
R  = diag([0.1, 0.1, 0.1, 1, 0.5]);

[H,f] = build_cost(F,G,dx0,Q,R,Qf,N);
[Aineq,bineq] = build_constraints(dx0, X_ref_win, U_ref_win);

options = optimoptions('quadprog','Display','off');

U = quadprog(H,f,Aineq,bineq,[],[],[],[],[],options);

u = U(1:nu);

end
