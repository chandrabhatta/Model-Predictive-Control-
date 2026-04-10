clear; clc; close all;

T = 200;
N = 30;
v_des = 1.0;
x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];

[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);
nx = size(xref, 1);
nu = size(uref, 1);

x0 = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

weight_sets = {
    diag([2000  2000  1000000 200 200 200  5  50]), diag([0.1, 0.1, 0.1, 10,  0.5]), 'Baseline';
    diag([10000 10000 1000000 200 200 200  5  50]), diag([0.1, 0.1, 0.1, 10,  0.5]), 'High Q_{pos}';
    diag([500   500   1000000 200 200 200  5  50]), diag([0.1, 0.1, 0.1, 10,  0.5]), 'Low Q_{pos}';
    diag([2000  2000  1000000 200 200 200 50  50]), diag([0.1, 0.1, 0.1, 10,  0.5]), 'High Q_{yaw}';
    diag([2000  2000  1000000 200 200 200  5  50]), diag([1.0, 1.0, 1.0, 50,  5.0]), 'High R';
};

colors = {'k', 'b', 'r', 'g', 'm'};

figure(1); hold on;
figure(2); hold on;

for idx = 1:size(weight_sets, 1)
    Q = weight_sets{idx, 1};
    R = weight_sets{idx, 2};
    lbl = weight_sets{idx, 3};

    [AdN, BdN] = linearize(xref(:,end), uref(:,end));
    Qf = dare(AdN, BdN, Q, R);

    x = x0;
    X_hist = zeros(nx, T);
    U_hist = zeros(nu, T);

    for k = 1:T
        i_end = min(k + N - 1, T);
        n_avail = i_end - k + 1;
        X_win = xref(:, k:i_end);
        U_win = uref(:, k:i_end);
        if n_avail < N
            X_win = [X_win, repmat(xref(:,end), 1, N - n_avail)];
            U_win = [U_win, repmat(uref(:,end), 1, N - n_avail)];
        end
        u = mpc_controller(x, xref(:,k), X_win, U_win, Qf, Q, R);
        x = rk4Integrator(x, u);
        X_hist(:,k) = x;
        U_hist(:,k) = u;
    end

    pos_err = sqrt((X_hist(1,:)-xref(1,:)).^2 + (X_hist(2,:)-xref(2,:)).^2);
    fprintf('%-18s  RMS = %.4f m\n', lbl, sqrt(mean(pos_err.^2)));

    figure(1);
    plot(tq, pos_err, colors{idx}, 'LineWidth', 1.5, 'DisplayName', lbl);

    figure(2);
    plot(X_hist(1,:), X_hist(2,:), colors{idx}, 'LineWidth', 1.5, 'DisplayName', lbl);
end

figure(1);
xlabel('Time [s]'); ylabel('Position error [m]');
title('Weight tuning: Position error'); legend show; grid on;

figure(2);
plot(xref(1,:), xref(2,:), 'k--', 'LineWidth', 1.2, 'DisplayName', 'Reference');
xlabel('X [m]'); ylabel('Y [m]');
title('Weight tuning: X and Y Trajectory'); legend show; grid on; axis equal;
