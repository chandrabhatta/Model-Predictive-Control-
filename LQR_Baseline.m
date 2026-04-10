%% --- Generate reference trajectory ---
N = 200;
v_des = 1.0;

x_wp = [0 2 4 6];
y_wp = [0 1 2 0];

[x_ref, u_ref, tq] = reference_generator(x_wp, y_wp, v_des, N);

nx = size(x_ref,1);
nu = size(u_ref,1);

%% --- Linearize along trajectory ---
Ad_seq = cell(1,N);
Bd_seq = cell(1,N);

for k = 1:N
    [Ad_seq{k}, Bd_seq{k}, ~] = linearize(x_ref(:,k), u_ref(:,k));
end

%% --- Define LQR weights ---
Q  = diag([2000 2000 1000000 200 200 200 5 50]);
R  = diag([0.1, 0.1, 0.1, 10, 0.5]);
[AdN, BdN] = linearize(x_ref(:,end), u_ref(:,end));
Qf = dare(AdN, BdN, Q, R);

%% --- Backward Riccati recursion ---
P = cell(1,N+1);
K = cell(1,N);

P{N+1} = Qf;

for k = N:-1:1
    Ad = Ad_seq{k};
    Bd = Bd_seq{k};
    
    K{k} = (R + Bd' * P{k+1} * Bd) \ (Bd' * P{k+1} * Ad);
    P{k} = Q + Ad' * P{k+1} * (Ad - Bd * K{k});
end

%% --- Forward simulate LQR on nonlinear dynamics ---
x = zeros(nx,N);
u = zeros(nu,N);

x(:,1) = x_ref(:,1) + [0.20; -0.10; 0; 0.05; 0; 0; 0.05; 0];

for k = 1:N-1
    u(:,k) = u_ref(:,k) - K{k} * (x(:,k) - x_ref(:,k));
    x(:,k+1) = rk4Integrator(x(:,k), u(:,k));
end

% last control
u(:,N) = u_ref(:,N) - K{N} * (x(:,N) - x_ref(:,N));

%% --- Plot reference vs LQR trajectory ---
figure;
plot3(x_ref(1,:), x_ref(2,:), x_ref(3,:), 'r--','LineWidth',2); hold on;
plot3(x(1,:), x(2,:), x(3,:), 'b','LineWidth',2); hold on;
%plot(u_ref(1,:), u_ref(3,:),'b','LineWidth',2); hold on;
%plot(u(1,:), u(3,:),'r','LineWidth',2);
xlabel('X'); ylabel('Y'); zlabel('Z');
%xlabel('U1'); ylabel('U3');
legend('Reference','LQR Tracking');
title('TV-LQR Tracking of Lunar Rover Trajectory');
grid on;