function [R,J] = rod(u,theta)
    %ROD Rodrigues formula: from axis-angle to rotation matrix
    
    %if theta > eps
    if true
        N =  [0 -u(3) u(2)
            u(3) 0 -u(1)
            -u(2) u(1) 0 ];
        R = eye(3) + N*sin(theta) + N^2*(1-cos(theta));
        J  = jacobianRotAx(u,theta);
    else
        R = eye(3); % degenerate case
    end
    
    
end


function J  = jacobianRotAx(u,th)
    %JACOBIANROTAX Jacobian of a rotation given by axis and angle
    
    % u asse per R
    % th angolo per R
    % out: matrice jacobiana 3x4 delle derivate
    %     (3 componenti di M x  4 parametri R asse+angolo)
    
    
    u1=u(1); u2=u(2); u3=u(3);
    
    
    N =  [0 -u3 u2;
        u3 0 -u1;
        -u2 u1 0;];
    
    % derivata risp. angolo
    J_th = N*( eye(3)*cos(th) + N*sin(th) ) ;
    
    % derivata risp. asse componente x
    A = [     0,    u2,    u3
        u2, -2*u1,     0
        u3,     0, -2*u1];
    B =[  0,  0,  0
        0,  0, -1
        0,  1,  0];
    J_u1  = ( sin(th)*B + (1-cos(th))*A ) ;
    
    % derivata risp. asse componente y
    A =[ -2*u2,    u1,     0
        u1,     0,    u3
        0,    u3, -2*u2];
    B =[  0,  0,  1
        0,  0,  0
        -1,  0,  0];
    J_u2  = (  sin(th)*B + (1-cos(th))*A ) ;
    
    % derivata risp. asse componente z
    A =[ -2*u3,     0,    u1
        0, -2*u3,    u2
        u1,    u2,     0];
    B =[  0, -1,  0
        1,  0,  0
        0,  0,  0];
    J_u3  = (  sin(th)*B + (1-cos(th))*A ) ;
    
    
    % metto tutto assieme
    J = [J_u1(:), J_u2(:), J_u3(:), J_th(:)];
    
    
end
