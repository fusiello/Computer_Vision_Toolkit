function [F12,e21] = fund(P1,P2)
%FUND Build fundamental matrix from camera matrices
% With 2 output parameters it also returns the epipole

Q = P2(:,1:3)/P1(:,1:3);
e21 = P2(:,4) - Q*P1(:,4);
F12=skew(e21)*Q;

% F12 is s.t. m2' * F12 * m1 = 0  
% the epipole belongs to image 2


