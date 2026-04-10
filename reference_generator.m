%% Generate a reference trajectory for MPC

function [xref, uref, tq] = reference_generator(x, y, v_des, N)
    % Defining waypoints:

    %x = [0 2 1.6 4.3];
    %y = [0 1 3.6 3.10];

    % Distance (m), time (s), Speed (m/s), angle (rad), 

    % --- Arc-length based timing ---
    dx = diff(x); % difference between waypoints
    dy = diff(y);
    dist = sqrt(dx.^2 + dy.^2); % distance between waypoints

    arc_length = [0 cumsum(dist)]; % cumulutively sum the distances
    %v_des = 2; % Constant speed of 2 m/s is assumed for the moon rover
    t = arc_length / v_des; % time steps

    % --- Spline construction ---
    ppx = spline(t, x);
    ppy = spline(t, y);

    % --- Query time ---
    tq = linspace(0, t(end), N);

    % --- Position ---
    xq = ppval(ppx, tq);
    yq = ppval(ppy, tq);
    zq = lunarTerrain(xq,yq);

    % --- Velocity ---
    ppx_dot = fnder(ppx);
    ppy_dot = fnder(ppy);

    vx = ppval(ppx_dot, tq);
    vy = ppval(ppy_dot, tq);

    % vertical component:

    vz = gradient(zq, tq);

    % --- Heading (Yaw) ---
    psi = atan2(vy, vx);
    psi = unwrap(psi);

    % --- Theta (Pitch) ---
    theta = atan2(vz, sqrt(vx.^2 + vy.^2));

    % Visualize the lunar terrain:

    [X, Y] = meshgrid( ...
        linspace(min(xq)-1, max(xq)+1, 100), ...
        linspace(min(yq)-1, max(yq)+1, 100));

    % Evaluate the terrain 
    Z = lunarTerrain(X,Y);

    % Input reference:

    % --- Accelerations ---
    ax = gradient(vx, tq);
    ay = gradient(vy, tq);
    az = gradient(vz, tq) + 1.62;   % add lunar gravity compensation (g_moon = 1.62 m/s^2)

    % Yaw rate:
    psi_dot = gradient(psi, tq);
    wbz = psi_dot;

    % Pitch rate:
    theta_dot = gradient(theta, tq);
    
    wby = theta_dot;

    figure;
    surf(X, Y, Z, 'EdgeColor', 'none');
    colormap(gray);
    hold on;

    plot3(xq, yq, zq, 'r', 'LineWidth', 2);

    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Lunar Terrain with Rover Trajectory');

    view(45, 30);
    axis tight;
    grid on;

    shading interp;
    
    camlight;
    lighting gouraud;

    % --- State reference ---
    xref = [xq;
            yq;
            zq;
            vx;
            vy;
            vz;
            psi;
            theta];

    % --- Input reference ---
    uref = [ax;
            ay;
            az;
            wby;
            wbz];

    % --- Add terrain map Pz = h(x,y) ---

end

