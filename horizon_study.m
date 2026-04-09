clear; clc; close all;

T = 200;
v_des = 1.0;
x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];

[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);
nx = size(xref, 1);

Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);

horizons = [5, 10, 20, 30, 40];
colors = {'r', 'b', 'g', 'm', 'k'};
rms_pos = zeros(1, length(horizons));
comp_time = zeros(1, length(horizons));

x0 = xref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

figure; hold on;
for idx = 1:length(horizons)
    N = horizons(idx);
    [AdN, BdN] = linearize(xref(:,end), uref(:,end));
    Qf = dare(AdN, BdN, Q, R);
    x = x0;
    X_hist = zeros(nx, T);
    t_start = tic;
    for k = 1:T
        i_end   = min(k + N - 1, T);
        n_avail = i_end - k + 1;
        X_win   = xref(:, k:i_end);
        U_win   = uref(:, k:i_end);
        if n_avail < N
            X_win = [X_win, repmat(xref(:,end), 1, N - n_avail)];
            U_win = [U_win, repmat(uref(:,end), 1, N - n_avail)];
        end
        u = mpc_controller(x, xref(:,k), X_win, U_win, Qf);
        x = rk4Integrator(x, u);
        X_hist(:,k) = x;
    end

    comp_time(idx) = toc(t_start);
    pos_err = sqrt((X_hist(1,:)-xref(1,:)).^2 + (X_hist(2,:)-xref(2,:)).^2);
    rms_pos(idx) = sqrt(mean(pos_err.^2));
    plot(tq, pos_err, colors{idx}, 'LineWidth', 1.5, 'DisplayName', sprintf('N = %d', N));
    fprintf('N = %2d  ,  RMS = %.4f m  ,  Time = %.2f s\n', N, rms_pos(idx), comp_time(idx));
end
xlabel('Time [s]'); ylabel('Position error [m]');
title('Horizon study: Position error'); legend show; grid on;

figure;
subplot(1,2,1);
bar(horizons, rms_pos, 'FaceColor', [0.2 0.5 0.8]);
xlabel('Prediction horizon N'); ylabel('RMS Position error [m]');
title('RMS error vs horizon'); grid on;

subplot(1,2,2);
bar(horizons, comp_time, 'FaceColor', [0.8 0.4 0.2]);
xlabel('Prediction horizon N'); ylabel('Computation time [s]');
title('Computation time vs horizon'); grid on;
