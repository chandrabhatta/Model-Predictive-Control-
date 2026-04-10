function [xref, uref, tq] = reference_generator(x, y, v_des, N)

% --- Arc-length timing ---
dx = diff(x);
dy = diff(y);
dist = sqrt(dx.^2 + dy.^2);

arc_length = [0 cumsum(dist)];
t = arc_length / v_des;

% --- spline trajectory ---
ppx = spline(t, x);
ppy = spline(t, y);

tq = linspace(0, t(end), N);

xq = ppval(ppx, tq);
yq = ppval(ppy, tq);
zq = lunarTerrain(xq, yq);

% --- velocities (from spline derivatives, smooth) ---
ppx_dot = fnder(ppx, 1);
ppy_dot = fnder(ppy, 1);

vx = ppval(ppx_dot, tq);
vy = ppval(ppy_dot, tq);

% smooth velocities (IMPORTANT FIX)
vx = smoothdata(vx, "sgolay", 7);
vy = smoothdata(vy, "sgolay", 7);

% --- vertical motion ---
vz = gradient(zq, tq);
vz = smoothdata(vz, "sgolay", 7);

% --- orientation ---
psi = atan2(vy, vx);
psi = unwrap(psi);

theta = atan2(vz, sqrt(vx.^2 + vy.^2));

% --- smoothed angular rates ---
psi_dot = gradient(psi, tq);
theta_dot = gradient(theta, tq);

psi_dot = smoothdata(psi_dot, "movmean", 5);
theta_dot = smoothdata(theta_dot, "movmean", 5);

% --- accelerations (smoothed instead of raw gradient spikes) ---
ax = gradient(vx, tq);
ay = gradient(vy, tq);

ax = smoothdata(ax, "movmean", 5);
ay = smoothdata(ay, "movmean", 5);

az = gradient(vz, tq) - 1.62;
az = smoothdata(az, "movmean", 5);

% --- inputs ---
uref = [ax;
        ay;
        az;
        theta_dot;
        psi_dot];

% --- IMPORTANT: clamp inputs to actuator limits ---
u_min = [-3; -2; -2; -1; -1];
u_max = [ 3;  2;  2;  1;  1];

uref = max(u_min, min(u_max, uref));

% --- states ---
xref = [xq;
        yq;
        zq;
        vx;
        vy;
        vz;
        psi;
        theta];

% --- plotting (unchanged) ---
figure;
[X, Y] = meshgrid( ...
    linspace(min(xq)-1, max(xq)+1, 100), ...
    linspace(min(yq)-1, max(yq)+1, 100));

Z = lunarTerrain(X,Y);

surf(X, Y, Z, 'EdgeColor', 'none');
colormap(gray); hold on;
plot3(xq, yq, zq, 'r', 'LineWidth', 2);

xlabel('X'); ylabel('Y'); zlabel('Z');
title('Lunar Terrain with Rover Trajectory');
view(45,30);
axis tight; grid on;
shading interp;
camlight; lighting gouraud;

end

