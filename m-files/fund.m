function [F12,e21] = fund(P1,P2)
%FUND Build fundamental matrix from camera matrices
% With 2 output parameters it also returns the epipole

F = zeros(3); e = zeros(3,1);

for i=1:3
    for j=1:3
        Q1=P1; Q1(j,:)=[];
        Q2=P2; Q2(i,:)=[]; 
        F(i,j)=(-1)^(i+j)*det([Q1;Q2]);
    end
    e(i) = det([P1(1:3,:); P2(i,:)]);
end

e21 = e./norm(e);
F12 = F./norm(F);


% Notation: F12 is s.t. m2' * F12 * m1 = 0  
% the epipole e21 belongs to image 2

% legacy method - assume invertible P1(:,1:3)
% Q = P2(:,1:3)/P1(:,1:3);
% e = P2(:,4) - Q*P1(:,4);
% F=skew(e)*Q;



