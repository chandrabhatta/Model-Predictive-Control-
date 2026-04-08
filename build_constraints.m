function [Aineq, bineq] = build_constraints(dx0, X_ref_win, U_ref_win)

%nx = size(X_ref_win,1); % calculate size of reference based on argument provided
nu = size(U_ref_win,1); % calculate size of input based on argument provided
N  = size(U_ref_win,2); % calculate the horizon length

% --- Input limits ---
u_min = [-2; -1; -0.5; -0.3; -0.5];
u_max = [ 2;  1;  0.5;  0.3;  0.5];

Umin = repmat(u_min, N, 1);
Umax = repmat(u_max, N, 1);

A_u = [eye(N*nu); -eye(N*nu)];
b_u = [Umax; -Umin];

% --- State limits ---
% Define reasonable limits (adjust as needed)
x_min = [-Inf; -Inf; -0.1; -2; -2; -0.5; -pi; -pi/6];
x_max = [ Inf;  Inf; 0.1;  2;  2;  0.5;  pi;  pi/6];

Xmin = repmat(x_min, N, 1);
Xmax = repmat(x_max, N, 1);

% --- Rollout prediction: F*dx0 + G*U gives predicted deviation ---
[F,G] = rolloutPrediction(X_ref_win, U_ref_win);

A_x = [ G; -G ];
b_x = [ Xmax - F*dx0; -Xmin + F*dx0 ];

Aineq = [A_u; A_x];
bineq = [b_u; b_x];

end
