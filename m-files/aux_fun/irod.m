function [u,theta] = irod(R)
    %IROD Inverse Rodrigues formula: from rotation matrix to axis-angle
    
    theta = real(acos((trace(R)-1)/2));
    u = [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
    
    % handle particular cases
    if (nnz(R-R')==0) % R-R'=0
        if (nnz(R-eye(3))==0) % R = I
            theta=0;
            u=[1;1;1]; % u can be any vector (the angle is zero!)
        else %R+I=2*u*u' has rank 1
            theta=pi;
            A=R+eye(3);
            u=A(:,1); % u can be computed by normalizing any column of R+I
        end
    end
    
    u = u/norm(u);
end
