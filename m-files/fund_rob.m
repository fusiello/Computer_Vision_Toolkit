function   [F, inliers] = fund_rob(m2,m1,varargin)
    %FUND_ROB  Robust fundamental matrix
    
    s = 8; % MSS size
    
    sampson_d = @(F,x)...
    sqrt(sum(reshape(sampson_fund(F,x(1:2,:),x(3:4,:)),4,[]).^2))';
    
    [~,inliers] = robust_fit([m1;m2],@(x,w)...
        fund_lin(x(3:4,:),x(1:2,:),w),sampson_d,s,varargin{:});
    
    if sum(inliers) < s
        warning('Fund_rob: not enough inliers');
    end
    
    % final linear fit
    F = fund_lin(m2(:,inliers), m1(:,inliers));
end

