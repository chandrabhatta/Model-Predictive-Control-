function [F,G] = build_prediction(Ad,Bd,N)

nx = size(Ad,1);
nu = size(Bd,2);

F = zeros(nx*N,nx);
G = zeros(nx*N,nu*N);

for i = 1:N
    F((i-1)*nx+1:i*nx,:) = Ad^i;
    
    for j = 1:i
        G((i-1)*nx+1:i*nx,(j-1)*nu+1:j*nu) = Ad^(i-j)*Bd;
    end
end

end