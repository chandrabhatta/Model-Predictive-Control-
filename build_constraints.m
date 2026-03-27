function [Aineq,bineq] = build_constraints(x0, xref, uref)

nx = 8;   % nx = 8
nu = 5;   % nu = 5
N  = 100; % N = 100

% Input constraints:
u_min = [-2; -1; -0.5; -0.3; -0.5];
u_max = [ 2;  1;  0.5;  0.3;  0.5];

Au = [eye(nu); -eye(nu)];
bu = [u_max; -u_min];

% Stack over the horizon:

Eu = kron(eye(N),Au);
bu_stack = repmat(bu, N, 1);

% State constraints:
x_min = [-1; -1; 0;   -1;  -0.5; -0.2; -0.26; -pi];
x_max = [ 5;  5; 1;    1;   0.5;  0.2;  0.26;  pi];

% Stack state constraints over the horizon:
Fx = kron(eye(N), [eye(nx); -eye(nx)]);
bx = repmat([x_max; -x_min], N, 1);

% Linearize crater obstacles:

crater_centers = [1 1.5; 3 2.5];
crater_radius  = [0.4, 0.3];

n_craters = size(crater_centers,1);

F_obs = [];
b_obs = [];

for k = 1:N

    px0 = xref(1);
    py0 = xref(2);
    
    for i = 1:n_craters
        
        xc = crater_centers(i,1);
        yc = crater_centers(i,2);
        r  = crater_radius(i);
        
        % Gradient (normal direction)
        a = px0 - xc;
        b = py0 - yc;
        
        % Linearized constraint
        c = a*px0 + b*py0 - r^2;
        
        % Build row for full state vector
        F_row = zeros(1, nx*N);
        
        idx = (k-1)*nx + 1;
        F_row(idx)   = a; % Px
        F_row(idx+1) = b; % Py
        
        F_obs = [F_obs; F_row];
        b_obs = [b_obs; c];
    end
end

F_total = [Fx; F_obs];
b_total = [bx; b_obs];


[P, S] = rolloutPrediction(xref, uref);

Aineq = [F_total * S;
         Eu];

bineq = [b_total - F_total * P * x0;
         bu_stack];

end
