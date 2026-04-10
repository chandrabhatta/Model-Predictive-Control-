function [K, alpha] = compute_terminal_set(Ad, Bd, Q, R, Qf, u_min, u_max)

K = dlqr(Ad, Bd, Q, R);
P = Qf;
P_inv = inv(P);

nu = size(K,1);
alpha = Inf;

for i = 1:nu
    K_i = K(i,:);

    delta_u = min(u_max(i), -u_min(i));

    if delta_u > 0
        alpha_i = delta_u^2 / (K_i * P_inv * K_i');
        alpha = min(alpha, alpha_i);
    end
end
end