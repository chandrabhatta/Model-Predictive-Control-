function u = lqr_controller(x, xr, ur, K)
dx = x - xr;
dx(7) = atan2(sin(dx(7)), cos(dx(7)));
u = ur - K * dx;
end
