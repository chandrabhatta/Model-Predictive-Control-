function [H, f] = build_cost(F, G, dx0, Q, R, Qf, U_ref_stack)


[~, Nnu]  = size(G);
N  = Nnu / size(R,1);

%These are introduced to add soft terrain constraints (REMOVE IF NEEDED):
rho = 1e4;  

% Build block-diagonal weighting matrices
if N == 1
    Qbar = Qf;
else
    Qbar = blkdiag(kron(eye(N-1), Q), Qf);
end
Rbar = kron(eye(N), R);

% Quadratic cost
H = G' * Qbar * G + Rbar;

% Linear term with reference tracking
f = G' * Qbar * (F*dx0);

% Soft terrain constraints:
H = blkdiag(H, rho * eye(N));
f = [f; zeros(N,1)];



% Ensure symmetry of H for quadprog
H = (H + H')/2;

end
