function [R,J] = eul(a)
    %EUL Rotation matrix from Euler angles
    
    % a contains the Euler anlges omega, phi, kappa
    
    R1 =  [1     0         0
        0 cos(a(1))  -sin(a(1))
        0 sin(a(1))  cos(a(1))]; % omega
    
    R2 = [cos(a(2))  0 sin(a(2))
        0          1       0
        -sin(a(2)) 0  cos(a(2))]; % phi
    
    R3 = [cos(a(3)) -sin(a(3))  0
        sin(a(3))   cos(a(3))  0
        0            0         1];  % kappa
    
    R = R3*R2*R1;
    
    if nargout > 1 % jacobian required
        
        % derivative wrt omega
        J_1 = R3 * R2* [ 0 0 0; 0 -sin(a(1)) -cos(a(1)); 0 cos(a(1))  -sin(a(1))];
        % derivative wrt phi
        J_2 = R3 * [-sin(a(2)) 0 cos(a(2)); 0 0 0; -cos(a(2)) 0 -sin(a(2)) ]  * R1 ;
        % derivative wrt kappa
        J_3 = [-sin(a(3)) -cos(a(3)) 0; cos(a(3)) -sin(a(3)) 0; 0 0 0] * R2 * R1;
        
        J = [J_1(:), J_2(:), J_3(:)] ;
    end
    
    if  abs(abs(a(2)) - pi/2) < 0.01 
        warning('Eulear rotation close to singularity')
    end
    
end
