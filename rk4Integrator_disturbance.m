function x_next = rk4Integrator_disturbance(x, u, d)

Ts = 0.05;

k1 = dynamics(x, u);
k2 = dynamics(x + (Ts/2)*k1, u);
k3 = dynamics(x + (Ts/2)*k2, u);
k4 = dynamics(x + Ts*k3, u);

x_next = x + (Ts/6) * (k1 + 2*k2 + 2*k3 + k4);

% Disturbance
if nargin < 3
    d = zeros(size(x));  
end
x_next = x_next + d;

% Angle wrapping
x_next(7) = atan2(sin(x_next(7)), cos(x_next(7)));
x_next(8) = atan2(sin(x_next(8)), cos(x_next(8)));

end