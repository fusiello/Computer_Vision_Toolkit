function R = rotation_synch(Z,A)
%ROTATION_SYNCH Rotation syncronization
    
    n = size(A,1);
    
    iD = diag(1./sum(A,2));             % inverse degree matrix
    [Q,~] = eigs( kron(iD,eye(3))*Z,3,'lr'); % top 3 eigenvectors
    
    Q = real(Q/(Q(1:3, 1:3))); % normalize first rotation to I
    % and guard against spurious complex valuess from roundoff
    
    R=cell(1,n);
    for i=1:n  % Projection onto SO(3)
        [U,~,V] = svd(Q(3*i-2:3*i,:));
        R{i} = U*diag([1,1,det(U*V')])*V';
    end
end

