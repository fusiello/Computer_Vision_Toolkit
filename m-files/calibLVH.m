function K  = calibLVH( H )
%CALIBLVH Camera calibration from infnity homographies
% Luong-Vieville-Hartley
      
    numV = length(H);
    S = duplication(3);
    
    L = [];                % build linear system
    for i =1:numV
        L = [L; (eye(9)-kron(H{i},H{i}))*S];
    end
    
    [~,~,V] = svd(L);       % solve
    B = inv(reshape(S*V(:,end),3,[]));
     
     if trace(B)<0 
        B=-B;% either B or -B is p.d.
    end  
    iK = diag([-1,-1,1])*chol(B); 
    
    K = inv(iK); K = K./K(3,3);
end

