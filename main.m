clear; clc; close all;

%% Simulation parameters
T  = 200;       % total time steps
Ts = 0.05;      % time step
N  = 30;        % prediction horizon

%% Waypoints for trajectory
x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

%% Generate reference trajectory
[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);
nx = size(xref,1);
nu = size(uref,1);

disp('Global reference input ranges:')

for i = 1:size(uref,1)
    fprintf('u%d: min = %.3f, max = %.3f\n', ...
        i, min(uref(i,:)), max(uref(i,:)));
end

%% Weights and terminal ingredients
Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);
[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);

u_min = [-3; -2; -2; -1; -1];
u_max = [ 3;  2;  2;  1;  1];
[K_terminal, alpha] = compute_terminal_set(AdN, BdN, Q, R, Qf, u_min, u_max, uref(:,end));

fprintf('Terminal Set \n');
fprintf('Alpha (terminal set size): %.6f\n', alpha);
Acl = AdN - BdN * K_terminal;
fprintf('Closed-loop eigenvalue magnitudes (should be < 1):\n');
disp(abs(eig(Acl))');

%% Initial state
x0 = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0]; % Ininital disturbance
x  = x0;

%% Storage
X_hist = zeros(nx, T);
U_hist = zeros(nu, T);

%% MPC simulation loop
for k = 1:T
    % Define current prediction window with end-of-trajectory padding
    i_end   = min(k + N - 1, T);
    n_avail = i_end - k + 1;
    X_win   = xref(:, k:i_end);
    U_win   = uref(:, k:i_end);
    if n_avail < N
        X_win = [X_win, repmat(xref(:,end), 1, N - n_avail)];
        U_win = [U_win, repmat(uref(:,end), 1, N - n_avail)];
    end

    xr = xref(:, k);

    % Compute MPC control input
    u = mpc_controller(x, xr, X_win, U_win, Qf);

    % Apply control input using RK4 integrator
    x = rk4Integrator(x, u, Ts);

    % Store results
    X_hist(:, k) = x;
    U_hist(:, k) = u;
end

%% Plot position
figure;
subplot(3,1,1);
plot(tq, X_hist(1,:), 'r', tq, xref(1,:), 'r--'); hold on;
plot(tq, X_hist(2,:), 'b', tq, xref(2,:), 'b--');
xlabel('Time [s]'); ylabel('Position [m]');
legend('X actual','X ref','Y actual','Y ref'); grid on;

%% Plot velocity
subplot(3,1,2);
plot(tq, X_hist(4,:), 'r', tq, xref(4,:), 'r--'); hold on;
plot(tq, X_hist(5,:), 'b', tq, xref(5,:), 'b--');
xlabel('Time [s]'); ylabel('Velocity [m/s]');
legend('Vx','Vx ref','Vy','Vy ref'); grid on;

%% Plot orientation
subplot(3,1,3);
plot(tq, X_hist(7,:), 'r', tq, xref(7,:), 'r--'); hold on;
plot(tq, X_hist(8,:), 'b', tq, xref(8,:), 'b--');
xlabel('Time [s]'); ylabel('Orientation [rad]');
legend('Psi','Psi ref','Theta','Theta ref'); grid on;

%% 3D trajectory
figure;
plot3(X_hist(1,:), X_hist(2,:), X_hist(3,:), 'b', 'LineWidth', 2); hold on;
plot3(xref(1,:), xref(2,:), xref(3,:), 'r--', 'LineWidth', 2);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
grid on; axis equal;
legend('Actual trajectory','Reference trajectory');
title('MPC trajectory tracking');

%% Performance summary
pos_err = sqrt((X_hist(1,:) - xref(1,:)).^2 + (X_hist(2,:) - xref(2,:)).^2);
yaw_err = abs(X_hist(7,:) - xref(7,:));
fprintf('Tracking performance \n');
fprintf('RMS position error: %.4f m\n', sqrt(mean(pos_err.^2)));
fprintf('Max position error: %.4f m\n', max(pos_err));
fprintf('RMS yaw error: %.4f deg\n', sqrt(mean(yaw_err.^2)) * 180/pi);

%% Terminal set membership check
dx_final = X_hist(:,end) - xref(:,end);
dx_final(7) = atan2(sin(dx_final(7)), cos(dx_final(7)));
V_final = dx_final' * Qf * dx_final;
fprintf('Terminal set check \n');
fprintf('V(x_T) = %.6f,  alpha = %.6f\n', V_final, alpha);
if V_final <= alpha
    fprintf('Final state is INSIDE the terminal set.\n');
else
    fprintf('Final state is OUTSIDE the terminal set.\n');
end

