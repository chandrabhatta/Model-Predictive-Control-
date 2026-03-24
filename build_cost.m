function [H,f] = build_cost(F,G,x0,Q,R,Qf,N)

nx = size(Q,1);
nu = size(R,1);

Qbar = kron(eye(N),Q);
Qbar(end-nx+1:end,end-nx+1:end) = Qf;

Rbar = kron(eye(N),R);

H = G'*Qbar*G + Rbar;
f = G'*Qbar*F*x0;

H = 2*H;
f = 2*f;

end