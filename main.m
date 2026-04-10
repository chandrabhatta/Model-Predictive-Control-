clear; clc; close all;

%% Simulation parameters
T  = 200;
Ts = 0.05;
N  = 15;

%% Waypoints
x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

%% Reference trajectory
[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);
nx = size(xref,1);
nu = size(uref,1);

%% Cost
Q = diag([20000 20000 1000000 200 200 200 5 50]);
R = diag([0.1, 0.1, 0.1, 10, 0.5]);

[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);

%% Initial state
x = xref(:,1);

%% Disturbance
d_true = [0; 0; 0; 0.05; -0.03; 0; 0.01; 0];

%% ===== EKF INITIALIZATION =====
nx_aug = 2*nx;

x_aug_hat = zeros(nx_aug,1);
P = 0.1 * eye(nx_aug);

Qk = 0.01 * eye(nx_aug);
Rk = 0.1 * eye(nx);

%% Storage
X_hist = zeros(nx, T);
U_hist = zeros(nu, T);
D_hist = zeros(nx, T);
D_true_hist = zeros(nx, T);

%% ===== MPC LOOP =====
for k = 1:T

    % Prediction window
    i = min(k+N-1, T);
    X_win = xref(:,k:i);
    U_win = uref(:,k:i);

    xr = xref(:,k);

    % Measurement
    y = x;

    % ===== EKF =====
    u_prev = U_hist(:, max(k-1,1));

    % Split
    xk = x_aug_hat(1:nx);
    dk = x_aug_hat(nx+1:end);

    % ---- Prediction ----
    x_pred = rk4Integrator(xk, u_prev, dk);
    d_pred = dk;

    x_aug_pred = [x_pred; d_pred];

    F = eye(nx_aug);
    F(1:nx, nx+1:end) = eye(nx);

    P_pred = F * P * F' + Qk;

    % ---- Update ----
    H = [eye(nx), zeros(nx)];

    K = P_pred * H' / (H * P_pred * H' + Rk);

    x_aug_hat = x_aug_pred + K*(y - H*x_aug_pred);

    P = (eye(nx_aug) - K*H) * P_pred;
    P = (P + P')/2;

    % Extract
    x_hat = x_aug_hat(1:nx);
    d_hat = x_aug_hat(nx+1:end);

    % ===== MPC =====
    u = mpc_controller(x_hat, xr, X_win, U_win, Qf);

    % System update
    x = rk4Integrator(x, u, d_true);

    % Store
    X_hist(:, k) = x;
    U_hist(:, k) = u;
    D_hist(:, k) = d_hat;
    D_true_hist(:, k) = d_true;
end

%% ===== PLOTS =====

figure;
subplot(3,1,1);
plot(tq, X_hist(1,:), 'r', tq, xref(1,:), 'r--'); hold on;
plot(tq, X_hist(2,:), 'b', tq, xref(2,:), 'b--');
legend('X','X ref','Y','Y ref'); grid on;

subplot(3,1,2);
plot(tq, X_hist(4,:), 'r', tq, xref(4,:), 'r--'); hold on;
plot(tq, X_hist(5,:), 'b', tq, xref(5,:), 'b--');
legend('Vx','Vx ref','Vy','Vy ref'); grid on;

subplot(3,1,3);
plot(tq, X_hist(7,:), 'r', tq, xref(7,:), 'r--'); hold on;
plot(tq, X_hist(8,:), 'b', tq, xref(8,:), 'b--');
legend('Psi','Psi ref','Theta','Theta ref'); grid on;

figure;
plot3(X_hist(1,:), X_hist(2,:), X_hist(3,:), 'b'); hold on;
plot3(xref(1,:), xref(2,:), xref(3,:), 'r--');
grid on; axis equal;
legend('Actual','Reference');
title('Trajectory');

figure;
subplot(3,1,1);
plot(tq, D_true_hist(4,:), 'k--'); hold on;
plot(tq, D_hist(4,:), 'b'); grid on;

subplot(3,1,2);
plot(tq, D_true_hist(5,:), 'k--'); hold on;
plot(tq, D_hist(5,:), 'b'); grid on;

subplot(3,1,3);
plot(tq, D_true_hist(7,:), 'k--'); hold on;
plot(tq, D_hist(7,:), 'b'); grid on;