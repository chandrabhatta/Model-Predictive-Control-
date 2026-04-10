clear; clc; close all;

disp('Running nominal MPC...');
[X_nom, U_nom, xref, tq] = run_mpc_nominal();

disp('Running EKF-based MPC...');
[X_ekf, U_ekf, D_ekf, D_true, ~, ~] = run_mpc_ekf();

nx = size(X_nom,1);

pos_err_nom = vecnorm(X_nom(1:2,:) - xref(1:2,:), 2, 1);
pos_err_ekf = vecnorm(X_ekf(1:2,:) - xref(1:2,:), 2, 1);

state_err_nom = vecnorm(X_nom - xref, 2, 1);
state_err_ekf = vecnorm(X_ekf - xref, 2, 1);

yaw_err_nom = abs(wrapToPi(X_nom(7,:) - xref(7,:)));
yaw_err_ekf = abs(wrapToPi(X_ekf(7,:) - xref(7,:)));

d_err = vecnorm(D_ekf - D_true, 2, 1);

fprintf('\n========== PERFORMANCE SUMMARY ==========\n');

fprintf('\n--- POSITION ERROR ---\n');
fprintf('Nominal MPC: RMSE = %.4f, MAX = %.4f\n', sqrt(mean(pos_err_nom.^2)), max(pos_err_nom));
fprintf('EKF MPC:     RMSE = %.4f, MAX = %.4f\n', sqrt(mean(pos_err_ekf.^2)), max(pos_err_ekf));

fprintf('\n--- FULL STATE ERROR ---\n');
fprintf('Nominal MPC: RMSE = %.4f\n', sqrt(mean(state_err_nom.^2)));
fprintf('EKF MPC:     RMSE = %.4f\n', sqrt(mean(state_err_ekf.^2)));

fprintf('\n--- YAW ERROR ---\n');
fprintf('Nominal MPC: RMSE = %.4f deg\n', sqrt(mean(yaw_err_nom.^2))*180/pi);
fprintf('EKF MPC:     RMSE = %.4f deg\n', sqrt(mean(yaw_err_ekf.^2))*180/pi);

fprintf('\n--- DISTURBANCE ERROR ---\n');
fprintf('EKF RMSE = %.5f\n', sqrt(mean(d_err.^2)));

figure;
subplot(2,1,1);
plot(tq, pos_err_nom, 'r', tq, pos_err_ekf, 'b', 'LineWidth', 1.5);
grid on;
legend('Nominal','EKF');
title('Position error');

subplot(2,1,2);
plot(tq, state_err_nom, 'r', tq, state_err_ekf, 'b', 'LineWidth', 1.5);
grid on;
legend('Nominal','EKF');
title('State error');

figure;
plot3(X_nom(1,:),X_nom(2,:),X_nom(3,:),'r', ...
      X_ekf(1,:),X_ekf(2,:),X_ekf(3,:),'b', ...
      xref(1,:),xref(2,:),xref(3,:),'k--');
grid on; axis equal;
legend('Nominal','EKF','Reference');
title('Trajectory comparison');

figure;
plot(tq,yaw_err_nom*180/pi,'r',tq,yaw_err_ekf*180/pi,'b','LineWidth',1.5);
grid on;
legend('Nominal','EKF');
title('Yaw error');

figure;
plot(tq,d_err,'k','LineWidth',1.5);
grid on;
title('Disturbance error');

function [X_hist, U_hist, xref, tq] = run_mpc_nominal()

T = 200; Ts = 0.05; N = 20;

x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

[xref, uref, tq] = reference_generator(x_waypoints, y_waypoints, v_des, T);

nx = size(xref,1);
nu = size(uref,1);

Q = diag([2000 2000 1e6 200 200 200 5 50]);
R = diag([0.1 0.1 0.1 10 0.5]);

[AdN,BdN] = linearize(xref(:,end), uref(:,end));
Qf = dare(AdN,BdN,Q,R);

K_terminal = compute_terminal_set(AdN,BdN,Q,R,Qf,[-4;-2;-2;-1;-1],[4;2;3;1;1],uref(:,end));

x = xref(:,1) + [0.2;-0.1;0;0.05;0;0;0.05;0];

G = [zeros(3,3); eye(3); zeros(2,3)];
%d_true = [0.05;-0.03;0.01];
d_true = 2 * [0.05; -0.03; 0.01];

X_hist = zeros(nx,T);
U_hist = zeros(nu,T);

for k = 1:T

    i_end = min(k+N-1,T);
    X_win = xref(:,k:i_end);
    U_win = uref(:,k:i_end);

    if size(X_win,2) < N
        X_win = [X_win repmat(xref(:,end),1,N-size(X_win,2))];
        U_win = [U_win repmat(uref(:,end),1,N-size(U_win,2))];
    end

    xr = xref(:,k);

    u = mpc_controller(x,xr,X_win,U_win,Qf,[],[],K_terminal);

    x = rk4Integrator(x,u) + G*d_true;

    X_hist(:,k) = x;
    U_hist(:,k) = u;
end

end

function [X_hist,U_hist,D_hist,D_true_hist,xref,tq] = run_mpc_ekf()

T = 200; Ts = 0.05; N = 30;

x_waypoints = [0 2 4 6];
y_waypoints = [0 1 2 0];
v_des = 1.0;

[xref,uref,tq] = reference_generator(x_waypoints,y_waypoints,v_des,T);

nx = size(xref,1);
nu = size(uref,1);

Q = diag([2000 2000 1e6 200 200 200 5 50]);
R = diag([0.1 0.1 0.1 10 0.5]);

[AdN,BdN] = linearize(xref(:,end),uref(:,end));
Qf = dare(AdN,BdN,Q,R);

K_terminal = compute_terminal_set(AdN,BdN,Q,R,Qf,[-4;-2;-2;-1;-1],[4;2;3;1;1],uref(:,end));

nd = 3;
nx_aug = nx + nd;

G = [zeros(3,3); eye(3); zeros(2,3)];

x = xref(:,1) + [0.2;-0.1;0;0.05;0;0;0.05;0];
d_true = [0.05;-0.03;0.01];

x_aug_hat = [x; zeros(nd,1)];
P = 0.1 * eye(nx_aug);

Qk = 0.01 * eye(nx_aug);
Rk = 0.1 * eye(nx);

X_hist = zeros(nx,T);
U_hist = zeros(nu,T);
D_hist = zeros(nd,T);
D_true_hist = repmat(d_true,1,T);

for k = 1:T

    i_end = min(k+N-1,T);
    X_win = xref(:,k:i_end);
    U_win = uref(:,k:i_end);

    if size(X_win,2) < N
        X_win = [X_win repmat(xref(:,end),1,N-size(X_win,2))];
        U_win = [U_win repmat(uref(:,end),1,N-size(U_win,2))];
    end

    xr = xref(:,k);
    y = x;

    if k == 1
        u_prev = U_win(:,1);
    else
        u_prev = U_hist(:,k-1);
    end

    xk = x_aug_hat(1:nx);
    dk = x_aug_hat(nx+1:end);

    [Ad,~,~] = linearize(xk,u_prev);

    x_nom = rk4Integrator(xk,u_prev);

    x_pred = x_nom + G*dk;
    d_pred = dk;

    x_aug_pred = [x_pred; d_pred];

    F = zeros(nx_aug,nx_aug);
    F(1:nx,1:nx) = Ad;
    F(1:nx,nx+1:end) = G;
    F(nx+1:end,nx+1:end) = eye(nd);

    P_pred = F*P*F' + Qk;

    H = [eye(nx) zeros(nx,nd)];

    S = H*P_pred*H' + Rk;
    K = P_pred*H'/S;

    x_aug_hat = x_aug_pred + K*(y - H*x_aug_pred);

    P = (eye(nx_aug)-K*H)*P_pred;

    x_hat = x_aug_hat(1:nx);
    d_hat = x_aug_hat(nx+1:end);

    u = mpc_controller(x_hat,xr,X_win,U_win,Qf,[],[],K_terminal);

    x = rk4Integrator(x,u) + G*d_true;

    X_hist(:,k) = x;
    U_hist(:,k) = u;
    D_hist(:,k) = d_hat;
end

end
