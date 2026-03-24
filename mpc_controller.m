function u = mpc_controller(x, xr, Ad, Bd, N)

dx0 = x - xr;

[F,G] = build_prediction(Ad,Bd,N);

nx = size(Ad,1);
nu = size(Bd,2);

Q  = eye(nx);
R  = 0.1*eye(nu);
Qf = Q;

[H,f] = build_cost(F,G,dx0,Q,R,Qf,N);
[Aineq,bineq] = build_constraints(nu,N);

options = optimoptions('quadprog','Display','off');

U = quadprog(H,f,Aineq,bineq,[],[],[],[],[],options);

u = U(1:nu);

end