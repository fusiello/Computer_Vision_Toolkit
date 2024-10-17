function [P,K] = calibSMZ( H )
%CALIBSMZ Camera calibration from object-image homographies
% Sturm-Maybank-Zhang

numV = length(H);
S = duplication(3);

L = [];                % build linear system
for i =1:numV
    L=[L;(kron(H{i}(:,2)',H{i}(:,1)'))*S;
        (kron(H{i}(:,1)',H{i}(:,1)')-kron(H{i}(:,2)',H{i}(:,2)'))*S];
end

[~,~,V] = svd(L);  % solve
B = reshape(S*V(:,end),3,[]);

if trace(B)<0
    B=-B; % either B or -B is p.d.
end
iK = diag([-1,-1,1])*chol(B); % force negative focal

K = inv(iK); K = K./K(3,3);

% compute extrinsic parameters
for i = 1:numV
    A = iK * H{i};
    A = A/norm(A(:,1)); % remove scale
    
    R = [A(:,1:2) , cross(A(:,1), A(:,2))];
    t = A(:,3);
    
    %  project onto SO(3)
    [U,~,V] = svd(R);
    R = U*diag([1,1,det(V'*U)])*V';
    
    if   [0,0,1] * (-R'*t)  < 0
        % change sign if cop is in the negative halfspace
        R(:,1:2) =  - R(:,1:2);
        t = - t;
    end
    P{i} =  K * [R,t];
end
end


