function   q = iquat(R)
    %IQUAT Computes the quaternion form of a rotation matrix
    
    % format: (scalar, vector)
    % BEWARE: the quaternion toolbox which uses (vector, scalar).
   
    
    % if abs(det(R) - 1.0)> 0.00001
    %     error('R must be a rotation')
    % end
    
    q(1)=0.5*sqrt(trace(R)+1);
    q(2)=0.5*sqrt(1+R(1,1)-R(2,2)-R(3,3));
    q(3)=0.5*sqrt(1-R(1,1)+R(2,2)-R(3,3));
    q(4)=0.5*sqrt(1-R(1,1)-R(2,2)+R(3,3));
    
    % to avoid division by 0
    [~,i]=max(q);
    switch i
        case 1
            q(2) = (R(3,2)-R(2,3))/(4*q(1));
            q(3) = (R(1,3)-R(3,1))/(4*q(1));
            q(4) = (R(2,1)-R(1,2))/(4*q(1));
        case 2
            q(1) = (R(3,2)-R(2,3))/(4*q(2));
            q(3) = (R(2,1)+R(1,2))/(4*q(2));
            q(4) = (R(1,3)+R(3,1))/(4*q(2));
        case 3
            q(1) = (R(1,3)-R(3,1))/(4*q(3));
            q(2) = (R(1,2)+R(2,1))/(4*q(3));
            q(4) = (R(2,3)+R(3,2))/(4*q(3));
        case 4
            q(1) = (R(2,1)-R(1,2))/(4*q(4));
            q(2) = (R(1,3)+R(3,1))/(4*q(4));
            q(3) = (R(2,3)+R(3,2))/(4*q(4));
    end
   
end
