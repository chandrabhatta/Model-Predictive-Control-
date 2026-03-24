clear
clc

nx = 8;
nu = 5;
N  = 10;
T  = 50;
dt = 0.1;

x = zeros(nx,1);
xr = zeros(nx,1);

Ad = eye(nx);
Bd = 0.1*randn(nx,nu);   % temporary placeholder

for k = 1:T
    
    u = mpc_controller(x,xr,Ad,Bd,N);
    
    % fake dynamics
    x = Ad*x + Bd*u;
    
    X(:,k) = x;
end

plot(X')
title('State evolution')