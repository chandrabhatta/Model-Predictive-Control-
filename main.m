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

u_min = [-4; -2; -2; -1; -1];
u_max = [ 4;  2;  3;  1;  1];  
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

    % Compute MPC control input (pass K_terminal for polytopic terminal set)
    u = mpc_controller(x, xr, X_win, U_win, Qf, [], [], K_terminal);

    % Apply control input using RK4 integrator
    x = rk4Integrator(x, u);

    % Store results
    X_hist(:, k) = x;
    U_hist(:, k) = u;
end

%% Plot position (X, Y, Z)
figure;
subplot(3,1,1);
plot(tq, X_hist(1,:), 'r',  tq, xref(1,:), 'r--'); hold on;
plot(tq, X_hist(2,:), 'b',  tq, xref(2,:), 'b--');
plot(tq, X_hist(3,:), 'g',  tq, xref(3,:), 'g--');
xlabel('Time [s]'); ylabel('Position [m]');
legend('X actual','X ref','Y actual','Y ref','Z actual','Z ref'); grid on;
title('Position tracking');

%% Plot velocity (Vx, Vy, Vz)
subplot(3,1,2);
plot(tq, X_hist(4,:), 'r',  tq, xref(4,:), 'r--'); hold on;
plot(tq, X_hist(5,:), 'b',  tq, xref(5,:), 'b--');
plot(tq, X_hist(6,:), 'g',  tq, xref(6,:), 'g--');
xlabel('Time [s]'); ylabel('Velocity [m/s]');
legend('Vx','Vx ref','Vy','Vy ref','Vz','Vz ref'); grid on;
title('Velocity tracking');

%% Plot orientation
subplot(3,1,3);
plot(tq, X_hist(7,:), 'r', tq, xref(7,:), 'r--'); hold on;
plot(tq, X_hist(8,:), 'b', tq, xref(8,:), 'b--');
xlabel('Time [s]'); ylabel('Orientation [rad]');
legend('Psi','Psi ref','Theta','Theta ref'); grid on;
title('Orientation tracking');

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

% Xf = { dx : u_min <= u_ref - K*dx <= u_max }
u_lqr_final = uref(:,end) - K_terminal * dx_final;
in_poly = all(u_lqr_final >= u_min) && all(u_lqr_final <= u_max);

fprintf('Terminal set check \n');
fprintf('LQR input at final state: [');
fprintf(' %.4f', u_lqr_final); fprintf(' ]\n');
fprintf('u_min:                    [');
fprintf(' %.4f', u_min); fprintf(' ]\n');
fprintf('u_max:                    [');
fprintf(' %.4f', u_max); fprintf(' ]\n');
if in_poly
    fprintf('Final state IS inside the terminal set (polytopic check passed).\n');
else
    fprintf('Final state is OUTSIDE the terminal set (polytopic check).\n');
    for i = 1:nu
        if u_lqr_final(i) < u_min(i)
            fprintf('  u%d = %.4f < u_min(%.2f)  [violation: %.4f]\n', ...
                i, u_lqr_final(i), u_min(i), u_min(i) - u_lqr_final(i));
        elseif u_lqr_final(i) > u_max(i)
            fprintf('  u%d = %.4f > u_max(%.2f)  [violation: %.4f]\n', ...
                i, u_lqr_final(i), u_max(i), u_lqr_final(i) - u_max(i));
        end
    end
end