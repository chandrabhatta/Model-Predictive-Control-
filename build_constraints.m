function [Aineq,bineq] = build_constraints(nu,N)

umin = -ones(nu,1);
umax =  ones(nu,1);

Umin = repmat(umin,N,1);
Umax = repmat(umax,N,1);

Aineq = [ eye(nu*N);
         -eye(nu*N)];

bineq = [Umax;
        -Umin];

end