clear; clc; close all;

T = 200;
N_mpc = 30;
v_des = 1.0;
x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];

[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);
nx = size(xref, 1);
nu = size(uref, 1);

Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);
[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);

x0 = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

%% MPC simulation
X_mpc = zeros(nx, T);
U_mpc = zeros(nu, T);
x = x0;
for k = 1:T
    i_end = min(k + N_mpc - 1, T);
    n_avail = i_end - k + 1;
    X_win = xref(:, k:i_end);
    U_win = uref(:, k:i_end);
    if n_avail < N_mpc
        X_win = [X_win, repmat(xref(:,end), 1, N_mpc - n_avail)];
        U_win = [U_win, repmat(uref(:,end), 1, N_mpc - n_avail)];
    end
    u = mpc_controller(x, xref(:,k), X_win, U_win, Qf);
    x = rk4Integrator(x, u);
    X_mpc(:,k) = x;
    U_mpc(:,k) = u;
end

%% TV-LQR simulation (backward Riccati then forward simulate)
Ad_seq = cell(1, T);
Bd_seq = cell(1, T);
for k = 1:T
    [Ad_seq{k}, Bd_seq{k}] = linearize(xref(:,k), uref(:,k));
end

P_lqr = cell(1, T+1);
K_lqr = cell(1, T);
P_lqr{T+1} = Qf;
for k = T:-1:1
    Ad = Ad_seq{k};  Bd = Bd_seq{k};
    K_lqr{k} = (R + Bd' * P_lqr{k+1} * Bd) \ (Bd' * P_lqr{k+1} * Ad);
    P_lqr{k} = Q + Ad' * P_lqr{k+1} * (Ad - Bd * K_lqr{k});
end

X_lqr = zeros(nx, T);
U_lqr = zeros(nu, T);
x = x0;
for k = 1:T
    dx = x - xref(:,k);
    dx(7) = atan2(sin(dx(7)), cos(dx(7)));
    u = uref(:,k) - K_lqr{k} * dx;
    x = rk4Integrator(x, u);
    X_lqr(:,k) = x;
    U_lqr(:,k) = u;
end

%% Results
err_mpc = sqrt((X_mpc(1,:)-xref(1,:)).^2 + (X_mpc(2,:)-xref(2,:)).^2);
err_lqr = sqrt((X_lqr(1,:)-xref(1,:)).^2 + (X_lqr(2,:)-xref(2,:)).^2);
fprintf('MPC  RMS position error: %.4f m\n', sqrt(mean(err_mpc.^2)));
fprintf('LQR  RMS position error: %.4f m\n', sqrt(mean(err_lqr.^2)));

figure;
subplot(2,1,1);
plot(tq, err_mpc, 'b', 'LineWidth', 1.5); hold on;
plot(tq, err_lqr, 'r--', 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('Position Error [m]');
legend('MPC', 'TV-LQR'); grid on;
title('Baseline Comparison: Position Error');

subplot(2,1,2);
plot(xref(1,:), xref(2,:), 'k--', 'LineWidth', 1.2); hold on;
plot(X_mpc(1,:), X_mpc(2,:), 'b', 'LineWidth', 1.5);
plot(X_lqr(1,:), X_lqr(2,:), 'r', 'LineWidth', 1.5);
xlabel('X [m]'); ylabel('Y [m]');
legend('Reference', 'MPC', 'TV-LQR'); grid on; axis equal;
title('Baseline Comparison: XY Trajectory');