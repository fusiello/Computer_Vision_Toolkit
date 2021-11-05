function [F,J] = par2F(f,U0,V0)
    %PAR2F Return a fundamental matrix and its jacobian from parameters
    
    % the parametrization is local around U0*V0'
    
    [U,JU] = eul(f(1:3));
    [V,JV] = eul(f(4:6));
    [s,Js] = ang2cart(f(7));
    D = diag([s;0]);
    % this 0-1 matrix is the jacobian of the diag([s;0]) operator
    JD = [ 1 0 0 0 0 0 0 0 0
           0 0 0 0 1 0 0 0 0]' *Js;
 
    F = U0*(U*D*V')*V0';
    
    K3 = commutation(3,3);
    J = kron(V0,U0)*[kron(V*D,eye(3))*JU, ...
        kron(eye(3),U*D)*K3*JV, kron(V,U)*JD ];
end





