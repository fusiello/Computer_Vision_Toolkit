function [H,J] = par2H(h)
    %PAR2H Return a homography matrix and its jacobian from parameters
    
    [v,Jv] = ang2cart(h(1:3));
    
    H(1:2,1:2) =  reshape(v,2,2);
    H(3,1) = h(4);
    H(3,2) = h(5);
    H(:,3) = h(6:8);
    
    % Jacobian 
    J = blkdiag(Jv,eye(5)); 
    % move row 5 in 3rd position
    J = J([1 2 5 3 4 6 7 8 9],:);
    
end
