function u = mpc_controller(x, xr, Ad, Bd, N)

dx0 = x - xr;

[F,G] = build_prediction(Ad,Bd,N);

nx = size(Ad,1);
nu = size(Bd,2);

Q  = diag([10, 10, 0.1, 1, 1, 0.1, 5, 1]);
R  = diag([0.1, 0.1, 0.1, 1, 0.5]);
Qf = dare(Ad, Bd, Q, R);

[H,f] = build_cost(F,G,dx0,Q,R,Qf,N);
[Aineq,bineq] = build_constraints(nu,N);

options = optimoptions('quadprog','Display','off');

U = quadprog(H,f,Aineq,bineq,[],[],[],[],[],options);

u = U(1:nu);

end
