function  R = quat(q)
     %QUAT Compute the rotation matrix given a unit quaternion
  
    % format: (scalar, vector)
    % BEWARE: the quaternion toolbox which uses (vector, scalar).
    
    
    % normalize the quaternion first
    q = q./norm(q);
    
    R = [ q(1)^2+q(2)^2-q(3)^2-q(4)^2,  2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3))
        2*(q(2)*q(3)+q(1)*q(4)), q(1)^2-q(2)^2+q(3)^2-q(4)^2,  2*(q(3)*q(4)-q(1)*q(2))
        2*(q(2)*q(4)-q(1)*q(3)),  2*(q(3)*q(4)+q(1)*q(2)),  q(1)^2-q(2)^2-q(3)^2+q(4)^2 ];
    
    
