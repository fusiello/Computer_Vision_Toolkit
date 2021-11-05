function K = fdh(H,u0,v0)
    %FDH compute focal distance from homography
    % The principal point (u0, v0) is needed
    
    S = duplication(3);
    
    L =[ (kron(H(:,2)',H(:,1)'))*S;
        (kron(H(:,1)',H(:,1)')-kron(H(:,2)',H(:,2)'))*S;
        0,   1, 0, 0,    0,  0;
        u0,  0, 1, 0,    0,  0;
        0,   0, 0, v0,   1,  0;
        u0^2,0, 0, v0^2, 0,  -1];
    
    % conditioning on rows and columns
    C = diag(1./sqrt(sum(L.^2,1)));
    R = diag(1./sqrt(sum(L.^2,2)));
    
    % soluzione del sistema lineare con precondizionamento
    b = C* ( (R * L * C ) \  R * [0;0;0;0;0;-1] ) ;
    
    B = reshape(S*b,3,[]);
    
    if trace(B)<0
        B=-B;% either B or -B is p.d.
    end
    iK = diag([-1,-1,1])*chol(B);
    
    K = inv(iK); K = K./K(3,3);
    
    %  Algorithm: Zhang 99 (MSR-TR-98-71)
end



