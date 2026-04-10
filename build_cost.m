function [H, f] = build_cost(F, G, dx0, Q, R, Qf, U_ref_stack)

% Inputs:
%   F, G          : rollout matrices
%   dx0           : current deviation x - x_ref
%   Q, R, Qf      : state/input weight matrices
%   U_ref_stack   : stacked reference input over horizon

[~, Nnu]  = size(G);
N  = Nnu / size(R,1);

%% These are introduced to add soft terrain constraints (REMOVE IF NEEDED):
rho = 1e4;   % penalty on slack (tune if needed)

%% Build block-diagonal weighting matrices
if N == 1
    Qbar = Qf;
else
    Qbar = blkdiag(kron(eye(N-1), Q), Qf);
end
Rbar = kron(eye(N), R);

% Quadratic cost
H = G' * Qbar * G + Rbar;

% Linear term with reference tracking
f = G' * Qbar * (F*dx0) - Rbar * U_ref_stack;

%% These are introduced to add soft terrain constraints (REMOVE IF NEEDED):
% --- Add slack variables (epsilon) ---
H = blkdiag(H, rho * eye(N));
f = [f; zeros(N,1)];



% Ensure symmetry of H for quadprog
H = (H + H')/2;

end
