function H = hom_synch(Z,A)
%HOM_SYNCH Homography syncronization
    
    n = size(A,1);
    
    iD = diag(1./sum(A,2));             % inverse degree matrix
    [Q,~]=eigs( kron(iD,eye(3))*Z,3,'lr'); % top 3 eigenvectors
    
    Q = real(Q/(Q(1:3, 1:3))); % normalize first homography to I
    % and guard against spurious complex valuess from roundoff
    
    H=cell(1,n);
    for i=1:n  % Force det = 1
        H{i} = Q(3*i-2:3*i,:)./nthroot(det(Q(3*i-2:3*i,:)),3);
    end
end


