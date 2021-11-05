function h = H2par(H)
    %H2PAR Homography matrix parametrization
    
    h = zeros(8,1);
    H = H./norm(H(1:2,1:2),'fro');
    
    h(1:3) = cart2ang(vec(H(1:2,1:2)));
    h(4) = H(3,1);
    h(5) = H(3,2);
    h(6:8) = H(:,3);
end

