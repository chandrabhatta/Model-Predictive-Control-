clear
clc

% We are making the robot follow a flat circular trajectory with a radius of 2 meters

nx = 8;
nu = 5;
N  = 15;
T  = 300;
Ts = 0.05;
R_traj = 2.0;
omega_c = 2*pi / (T * Ts);
g = 9.81;

X_ref = zeros(nx,T+1);
U_ref = zeros(nu,T);

for k = 1:T+1
    t_k = (k-1) * Ts;
    X_ref(1,k) = R_traj * cos(omega_c * t_k); % X
    X_ref(2,k) = R_traj * sin(omega_c * t_k); % Y
    X_ref(3,k) = 0; % Z
    X_ref(4,k) = -R_traj * omega_c * sin(omega_c * t_k); % Vx
    X_ref(5,k) = R_traj * omega_c * cos(omega_c * t_k); % Vy
    X_ref(6,k) = 0; % Vz
    X_ref(7,k) = atan2(X_ref(5,k), X_ref(4,k)); % psi
    X_ref(8,k) = 0; % theta
end
for k= 1:T
    U_ref(1,k) = 0;  % ax
    U_ref(2,k) = R_traj * omega_c^2; % ay
    U_ref(3,k) = g; % az
    U_ref(4,k) = 0; % omega_by
    U_ref(5,k) = omega_c; % omega_bz
end

Q = diag([10, 10, 0.1, 1, 1, 0.1, 5, 1]);
R = diag([0.1, 0.1, 0.1, 1.0, 0.5]);

[Ad0, Bd0] = linearize(X_ref(:,1), U_ref(:,1));
Qf = dare(Ad0, Bd0, Q, R);

x = X_ref(:,1) + [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]; % Put initial disturbance here
X_hist = zeros(nx, T);

for k = 1:T
    i = min(k + N, T + 1);
    n = i - k;
    if n >= N
        X_win = X_ref(:, k:k+N);
        U_win = U_ref(:, k:k+N-1);
    else
        X_win = [X_ref(:,k:i), repmat(X_ref(:,end), 1, N - n)]; % Filling the rest of the matrix with the same vector to avoid a crash
        U_win = [U_ref(:,k:min(k+N-1, T)), repmat(U_ref(:,end), 1, N - n)]; % Filling the rest of the matrix with the same vector to avoid a crash
    end

    xr = X_ref(:, k);
    u = mpc_controller(x, xr, X_win, U_win, Qf);
    x = rk4Integrator(x, u);
    X_hist(:, k) = x;
end
