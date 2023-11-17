function [R,t] = relative_lin(q2, q1, K2, K1)
%RELATIVE_LIN Relative orientation

q1 = K1\ensure_homogeneous(q1);
q2 = K2\ensure_homogeneous(q2);
% now q1 and q2 are homogeneous NIC

E = eight_points(q2(1:2,:), q1(1:2,:)) ;

[U,~,V] = svd(E);

S1 = [0 1 0; -1 0 0; 0 0 0];
R1 = [0 -1 0; 1 0 0; 0 0 1];

for j = 1:4
    % 4 combinations
    S = (-1)^j * U*S1*U';
    
    if j==3   R1 = R1'; end
    
    t=[S(3,2) S(1,3) S(2,1)]';
    R = det(U*V')*U*R1*V';
    
    % solve t = z2*q2 - z1*R*q1 for z1 and z2
    z = [ q2(:,1), -R*q1(:,1)] \t;
    % 3D points must be in front of both cameras
    if z(1)>0 && z(2)>0  break,  end
    
end


