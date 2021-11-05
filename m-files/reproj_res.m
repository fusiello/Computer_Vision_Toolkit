function [res,A,B,C,D]  = reproj_res(K,kappa,G,M,m)
    %REPROJ_RES Reprojection residual and jacobian for 1 point - 1 camera
    % There are two residual (u and v)
    % G = [R,t]
    
    M = [M(:);1];   % make homogeneous
    q = G * M;      
    
    % perspective division occurs here
    De = [ 1/q(3), 0, -q(1)/q(3)^2; 0, 1/q(3), -q(2)/q(3)^2];
    [qp,Dr,Dkappa] = rdx(kappa,q(1:2)./q(3),eye(3));
    
    res =  htx(K,qp) - m(:);
    
    A = K(1:2,1:2) * Dr * De * kron( M',eye(3) ); % wrt external par.
    B = K(1:2,1:2) * Dr * De * G(:,1:3);          % wrt 3D point
    C = kron([qp;1]', eye(3)); C(3,:)=[];         % wrt internal par.
    D = K(1:2,1:2) * Dkappa;                      % wrt distortion par.
end




%     Dfp = [ 1/W(3), 0, -W(1)/W(3)^2; 0, 1/W(3), -W(2)/W(3)^2; 0 0 0];
%     Dm(3,3) = 0; % pad Dm with zeros (homogeneous)
%
%     res =  htx(K,Mp) - m(:);
%
%     A = K * Dm * Dfp * kron( M',eye(3) ); A(3,:)=[];  % wrt external par.
%     B = K * Dm * Dfp * G(:,1:3);          B(3,:)=[];  % wrt 3D point
%     C = kron([Mp;1]', eye(3));            C(3,:)=[];  % wrt internal par.
%     D = K(1:2,1:2) * Dk;