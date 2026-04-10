clear; clc; close all;

%% Simulation parameters
T  = 200;
Ts = 0.05;
N  = 30;

%% Waypoints
x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

%% Reference trajectory
[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);
nx = size(xref,1);
nu = size(uref,1);

%% Cost
Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);

[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);

%% Terminal set 
u_min = [-4; -2; -2; -1; -1];
u_max = [ 4;  2;  3;  1;  1];

[K_terminal, alpha] = compute_terminal_set(AdN, BdN, Q, R, Qf, u_min, u_max, uref(:,end));

fprintf('Alpha: %.6f\n', alpha);

%% Initial state
x = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

%% True disturbance 
d_true = [0; 0; 0; 0.05; -0.03; 0; 0.01; 0];

%% Initialize EKF
nx_aug = 2*nx;

x_aug_hat = [x; zeros(nx,1)];
P = 0.1 * eye(nx_aug);

Qk = 0.01 * eye(nx_aug);
Rk = 0.1 * eye(nx);

%% Storage
X_hist = zeros(nx, T);
U_hist = zeros(nu, T);
D_hist = zeros(nx, T);
D_true_hist = zeros(nx, T);

%% MPC Loop
for k = 1:T

    % Prediction window (with padding)
    i_end   = min(k + N - 1, T);
    n_avail = i_end - k + 1;

    X_win = xref(:, k:i_end);
    U_win = uref(:, k:i_end);

    if n_avail < N
        X_win = [X_win, repmat(xref(:,end), 1, N - n_avail)];
        U_win = [U_win, repmat(uref(:,end), 1, N - n_avail)];
    end

    xr = xref(:, k);

    %% EKF
    y = x;

    if k == 1
        u_prev = U_win(:,1);
    else
        u_prev = U_hist(:,k-1);
    end

    % Split state
    xk = x_aug_hat(1:nx);
    dk = x_aug_hat(nx+1:end);

    % Prediction 
    x_pred = rk4Integrator(xk, u_prev, dk);
    d_pred = dk;

    x_aug_pred = [x_pred; d_pred];

    % Simple Jacobian
    F = eye(nx_aug);
    F(1:nx, nx+1:end) = eye(nx);

    P_pred = F * P * F' + Qk;

    % Update 
    H = [eye(nx), zeros(nx)];

    S = H * P_pred * H' + Rk;
    K = P_pred * H' / S;

    x_aug_hat = x_aug_pred + K*(y - H*x_aug_pred);

    P = (eye(nx_aug) - K*H) * P_pred;
    P = (P + P')/2;

    % Extract estimates
    x_hat = x_aug_hat(1:nx);
    d_hat = x_aug_hat(nx+1:end);

    %% MPC Using Estimate
    u = mpc_controller(x_hat, xr, X_win, U_win, Qf, [], [], K_terminal);

    x = rk4Integrator(x, u, d_true);

    %% Store
    X_hist(:, k) = x;
    U_hist(:, k) = u;
    D_hist(:, k) = d_hat;
    D_true_hist(:, k) = d_true;
end

%% Plots:

% Position
figure;
subplot(3,1,1);
plot(tq, X_hist(1,:), 'r', tq, xref(1,:), 'r--'); hold on;
plot(tq, X_hist(2,:), 'b', tq, xref(2,:), 'b--');
plot(tq, X_hist(3,:), 'g', tq, xref(3,:), 'g--');
legend('X','X ref','Y','Y ref','Z','Z ref'); grid on;
title('Position');

% Velocity
subplot(3,1,2);
plot(tq, X_hist(4,:), 'r', tq, xref(4,:), 'r--'); hold on;
plot(tq, X_hist(5,:), 'b', tq, xref(5,:), 'b--');
plot(tq, X_hist(6,:), 'g', tq, xref(6,:), 'g--');
legend('Vx','Vx ref','Vy','Vy ref','Vz','Vz ref'); grid on;
title('Velocity');

% Orientation
subplot(3,1,3);
plot(tq, X_hist(7,:), 'r', tq, xref(7,:), 'r--'); hold on;
plot(tq, X_hist(8,:), 'b', tq, xref(8,:), 'b--');
legend('Psi','Psi ref','Theta','Theta ref'); grid on;
title('Orientation');

% 3D trajectory
figure;
plot3(X_hist(1,:), X_hist(2,:), X_hist(3,:), 'b', 'LineWidth', 2); hold on;
plot3(xref(1,:), xref(2,:), xref(3,:), 'r--', 'LineWidth', 2);
grid on; axis equal;
legend('Actual','Reference');
title('Trajectory');
