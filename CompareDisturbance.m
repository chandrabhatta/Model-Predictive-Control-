clear; clc; close all;


disp('Running nominal MPC...');
[X_nom, U_nom, xref, tq] = run_mpc_nominal();

disp('Running EKF-based MPC...');
[X_ekf, U_ekf, D_ekf, D_true, ~, ~] = run_mpc_ekf();

nx = size(X_nom,1);

% Error metrics

% Position error (XY)
pos_err_nom = vecnorm(X_nom(1:2,:) - xref(1:2,:), 2, 1);
pos_err_ekf = vecnorm(X_ekf(1:2,:) - xref(1:2,:), 2, 1);

% Full-state error (L2)
state_err_nom = vecnorm(X_nom - xref, 2, 1);
state_err_ekf = vecnorm(X_ekf - xref, 2, 1);

% Yaw error (index 7)
yaw_err_nom = abs(wrapToPi(X_nom(7,:) - xref(7,:)));
yaw_err_ekf = abs(wrapToPi(X_ekf(7,:) - xref(7,:)));

% ===================== DISTURBANCE ESTIMATION QUALITY =====================

if exist('D_ekf','var')
    d_err = vecnorm(D_ekf - D_true, 2, 1);
else
    d_err = [];
end

% ===================== NUMERICAL SUMMARY =====================

fprintf('\n========== PERFORMANCE SUMMARY ==========\n');

fprintf('\n--- POSITION ERROR ---\n');
fprintf('Nominal MPC: RMSE = %.4f, MAX = %.4f\n', ...
    sqrt(mean(pos_err_nom.^2)), max(pos_err_nom));

fprintf('EKF MPC:     RMSE = %.4f, MAX = %.4f\n', ...
    sqrt(mean(pos_err_ekf.^2)), max(pos_err_ekf));

fprintf('\n--- FULL STATE ERROR ---\n');
fprintf('Nominal MPC: RMSE = %.4f\n', sqrt(mean(state_err_nom.^2)));
fprintf('EKF MPC:     RMSE = %.4f\n', sqrt(mean(state_err_ekf.^2)));

fprintf('\n--- YAW ERROR ---\n');
fprintf('Nominal MPC: RMSE = %.4f deg\n', sqrt(mean(yaw_err_nom.^2))*180/pi);
fprintf('EKF MPC:     RMSE = %.4f deg\n', sqrt(mean(yaw_err_ekf.^2))*180/pi);

if ~isempty(d_err)
    fprintf('\n--- DISTURBANCE ESTIMATION ---\n');
    fprintf('EKF disturbance RMSE = %.5f\n', sqrt(mean(d_err.^2)));
end

% Plots:

% 1. Position tracking comparison
figure;
subplot(2,1,1);
plot(tq, pos_err_nom, 'r', 'LineWidth', 1.5); hold on;
plot(tq, pos_err_ekf, 'b', 'LineWidth', 1.5);
grid on;
legend('Nominal MPC','EKF-MPC');
title('Position error (XY norm)');
ylabel('Error [m]');

subplot(2,1,2);
plot(tq, state_err_nom, 'r', 'LineWidth', 1.5); hold on;
plot(tq, state_err_ekf, 'b', 'LineWidth', 1.5);
grid on;
legend('Nominal MPC','EKF-MPC');
title('Full state error');
ylabel('||x - x_{ref}||');
xlabel('Time [s]');

% 2. 3D trajectory comparison
figure;
plot3(X_nom(1,:), X_nom(2,:), X_nom(3,:), 'r', 'LineWidth', 2); hold on;
plot3(X_ekf(1,:), X_ekf(2,:), X_ekf(3,:), 'b', 'LineWidth', 2);
plot3(xref(1,:), xref(2,:), xref(3,:), 'k--', 'LineWidth', 1.5);

grid on; axis equal;
legend('Nominal MPC','EKF MPC','Reference');
title('Trajectory comparison');
xlabel('X'); ylabel('Y'); zlabel('Z');

% 3. Yaw tracking comparison
figure;
plot(tq, yaw_err_nom*180/pi, 'r', 'LineWidth', 1.5); hold on;
plot(tq, yaw_err_ekf*180/pi, 'b', 'LineWidth', 1.5);
grid on;
legend('Nominal MPC','EKF MPC');
title('Yaw error');
ylabel('deg');
xlabel('Time [s]');

% 4. Disturbance estimation (if available)
if ~isempty(d_err)
    figure;
    plot(tq, d_err, 'k', 'LineWidth', 1.5);
    grid on;
    title('Disturbance estimation error (EKF)');
    xlabel('Time [s]');
    ylabel('||d_{hat} - d||');
end

function [X_hist, U_hist, xref, tq] = run_mpc_nominal()

clearvars -except X_hist U_hist xref tq   % optional if you reuse workspace

% Parematers:
T  = 200;
Ts = 0.05;
N  = 30;

x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

% Reference:
[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);

nx = size(xref,1);
nu = size(uref,1);

% Cost:
Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);

[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);

u_min = [-4; -2; -2; -1; -1];
u_max = [ 4;  2;  3;  1;  1];

[K_terminal, ~] = compute_terminal_set(AdN, BdN, Q, R, Qf, u_min, u_max, uref(:,end));

% Initial State:
x = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

% Storage:
X_hist = zeros(nx, T);
U_hist = zeros(nu, T);

% MPC Loop:
for k = 1:T

    % reference window:
    i_end = min(k + N - 1, T);
    n_avail = i_end - k + 1;

    X_win = xref(:, k:i_end);
    U_win = uref(:, k:i_end);

    if n_avail < N
        X_win = [X_win, repmat(xref(:,end), 1, N - n_avail)];
        U_win = [U_win, repmat(uref(:,end), 1, N - n_avail)];
    end

    xr = xref(:, k);

    % MPC control:
    u = mpc_controller(x, xr, X_win, U_win, Qf, [], [], K_terminal);

    x = rk4Integrator_disturbance(x, u, Ts);

    X_hist(:,k) = x;
    U_hist(:,k) = u;
end

end

function [X_hist, U_hist, D_hist, D_true_hist, xref, tq] = run_mpc_ekf()

% Parameters:
T  = 200;
Ts = 0.05;
N  = 30;

x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

% Reference:
[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);

nx = size(xref,1);
nu = size(uref,1);

% Cost:
Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);

[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);

u_min = [-4; -2; -2; -1; -1];
u_max = [ 4;  2;  3;  1;  1];

[K_terminal, ~] = compute_terminal_set(AdN, BdN, Q, R, Qf, u_min, u_max, uref(:,end));

% True System:
x = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

d_true = [0; 0; 0; 0.05; -0.03; 0; 0.01; 0];

% EKF
nx_aug = 2*nx;
x_aug_hat = [x; zeros(nx,1)];
P = 0.1 * eye(nx_aug);

Qk = 0.01 * eye(nx_aug);
Rk = 0.1 * eye(nx);

% Storage
X_hist = zeros(nx, T);
U_hist = zeros(nu, T);
D_hist = zeros(nx, T);
D_true_hist = zeros(nx, T);

% MPC Loop:
for k = 1:T

    % reference window:
    i_end = min(k + N - 1, T);
    X_win = xref(:, k:i_end);
    U_win = uref(:, k:i_end);

    if size(X_win,2) < N
        X_win = [X_win, repmat(xref(:,end), 1, N-size(X_win,2))];
        U_win = [U_win, repmat(uref(:,end), 1, N-size(U_win,2))];
    end

    xr = xref(:, k);

    % EKF
    y = x;

    if k == 1
        u_prev = U_win(:,1);
    else
        u_prev = U_hist(:,k-1);
    end

    xk = x_aug_hat(1:nx);
    dk = x_aug_hat(nx+1:end);

    % prediction:

    x_pred = rk4Integrator(xk, u_prev, dk);
    d_pred = dk;

    x_aug_pred = [x_pred; d_pred];

    F = eye(nx_aug);
    F(1:nx, nx+1:end) = eye(nx);

    P_pred = F * P * F' + Qk;

    % ---- update ----
    H = [eye(nx), zeros(nx)];
    S = H*P_pred*H' + Rk;
    K = P_pred*H'/S;

    x_aug_hat = x_aug_pred + K*(y - H*x_aug_pred);
    P = (eye(nx_aug)-K*H)*P_pred;

    x_hat = x_aug_hat(1:nx);
    d_hat = x_aug_hat(nx+1:end);

    % MPC:
    u = mpc_controller(x_hat, xr, X_win, U_win, Qf, [], [], K_terminal);

    x = rk4Integrator_disturbance(x, u, d_true);

    % Store
    X_hist(:,k) = x;
    U_hist(:,k) = u;
    D_hist(:,k) = d_hat;
    D_true_hist(:,k) = d_true;

end

end

