function   [H, inliers] = hom_rob(m2,m1,varargin)
    %HOM_ROB  Robust homography fit
    
    s = 4; % MSS size
    
    sampson_d = @(H,x)sqrt(sum(reshape(sampson_hom(H,x(1:2,:),x(3:4,:)),4,[]).^2))';
    
    [~,inliers] = robust_fit([m1;m2],@(x,w)hom_lin(x(3:4,:),x(1:2,:),w),...
        sampson_d, s, varargin{:});
    
    if sum(inliers) < s
        warning('Hom_rob: not enough inliers'); 
    end
    
    % final linear fit
    H = hom_lin(m2(:,inliers), m1(:,inliers));
end

% IRLS not implemented yet (weighted hom_lin would be needed)
