function [Phi, Gamma] = rolloutPrediction(X_ref_win, U_ref_win)

N  = 15;
nx = size(X_ref_win, 1);
nu = size(U_ref_win, 1);

Phi   = zeros(N*nx, nx);
Gamma = zeros(N*nx, N*nu);

Ad_seq = cell(N, 1);
Bd_seq = cell(N, 1);
for i = 1:N
    [Ad_seq{i}, Bd_seq{i}] = linearize(X_ref_win(:, i), U_ref_win(:, i));
end

A_prod = eye(nx);
for i = 1:N
    A_prod = Ad_seq{i} * A_prod;
    Phi((i-1)*nx+1 : i*nx, :) = A_prod;
end

for j = 1:N
    col_block = Bd_seq{j};
    for i = j:N
        Gamma((i-1)*nx+1 : i*nx, (j-1)*nu+1 : j*nu) = col_block;
        if i < N
            col_block = Ad_seq{i+1} * col_block;
        end
    end
end
end
