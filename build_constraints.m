function [Aineq, bineq] = build_constraints(dx0, X_ref_win, U_ref_win)

nu = size(U_ref_win,1); 
N  = size(U_ref_win,2); 

u_min = [-3; -2; -2; -1; -1];
u_max = [ 3;  2;  3;  1;  1];

Umin = repmat(u_min, N, 1);
Umax = repmat(u_max, N, 1);
U_ref_stack = reshape(U_ref_win, [], 1);

A_u = [eye(N*nu); -eye(N*nu)];
b_u = [Umax - U_ref_stack; -Umin + U_ref_stack];

% State limits 
x_min = [-Inf; -Inf; -Inf; -2; -2; -0.5; -pi; -pi/6];
x_max = [ Inf;  Inf; Inf;  2;  2;  0.5;  pi;  pi/6];


Xmin = repmat(x_min, N, 1);
Xmax = repmat(x_max, N, 1);

[F,G] = rolloutPrediction(X_ref_win, U_ref_win);

% These can be removed if the slack constraint doesn't work:

nx = size(X_ref_win,1);
z_idx = 3:nx:(N*nx);  

Gz = G(z_idx, :);
Fz = F(z_idx, :);

z_min = lunarTerrain(X_ref_win(1,:), X_ref_win(2,:));
z_min = reshape(z_min, [], 1);

margin = 0.05;
z_min = z_min - margin;

A_z = [-Gz, -eye(N)];
b_z = -z_min + Fz * dx0;

A_x = [ G; -G ];
b_x = [ Xmax - F*dx0; -Xmin + F*dx0 ];

Aineq = [A_u; A_x];
bineq = [b_u; b_x];

Aineq = [Aineq, zeros(size(Aineq,1), N)];
Aineq = [Aineq; A_z];
bineq = [bineq; b_z];

end
