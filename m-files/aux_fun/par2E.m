function [E,J] = par2E(g)
    %PAR2E Return an essential matrix and its jacobian from parameters
    
    % build E from parameters (5 angles)
    [R,DR] = eul(g(1:3));
    [t,Dt] = ang2cart(g(4:5));    
    [S,DS] = skew(t);
    E = S*R;

    J = [kron(eye(3),S) * DR , kron(R',eye(3)) * DS * Dt];
end
