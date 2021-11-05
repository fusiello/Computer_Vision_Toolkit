function H = hom_lin(m2,m1,w)
    %HOM_LIN  Homography with DLT algorithm
    
    if nargin < 3 || isempty(w)
        w = ones(size(m1,2),1); % weights
    end 
    
    % pre-conditioning
    [T1,m1] = precond(m1); [T2,m2] = precond(m2);
    
    H = dlt(m2, m1, w);
    
    % post-conditioning
    H = T2\H * T1;
end













