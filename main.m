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

%% Initial linearization to get terminal cost
%Q  = diag([10, 10, 1e-10, 1, 1, 0.1, 5, 1]); % 0.1, 0.1 at 6 originally
Q = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);
[AdN, BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN, BdN, Q, R);  % terminal cost now uses last point


%% Initial state
x0 = xref(:,1);  % can add initial disturbance
x  = x0;

%% Storage
X_hist = zeros(nx, T);
U_hist = zeros(nu, T);

%% MPC simulation loop
for k = 1:T
    % Define current prediction window
    i = min(k+N-1, T);
    X_win = xref(:, k:i);
    U_win = uref(:, k:i);

    % Current reference state
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
title('MPC Trajectory Tracking');

