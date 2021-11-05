function [F,e2] = fund(P1,P2)
%FUND Build fundamental matrix from camera matrices
% With 2 output parameters it also returns the epipole

Q = P2(:,1:3)/P1(:,1:3);
e2 = P2(:,4) - Q*P1(:,4);

F=skew(e2)*Q;



