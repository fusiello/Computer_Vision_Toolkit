function [H1,H2, K] = rectifyF(m1, m2, sz, f0)
    % RECTIFYF Epipolar rectification (uncalibrated)
  
    m1=ensure_homogeneous(m1);
    m2=ensure_homogeneous(m2);
    
    % first optimize with fixed f (to avoid jacobian singularity)
    p_out = lsq_nonlin(@(x)costRectif(x,sz,m1,m2), [0 0 0 0 0]);
      
    % then add f to the optimization
    p_out = lsq_nonlin(@(x)costRectif(x,sz,m1,m2), ...
        [p_out(:); log(3)*log(f0/(sz(1)+sz(2)))]);
    
    % recover K
    f = 3^p_out(6)*(sz(1)+sz(2));
    K = [f,0,sz(1)/2; 0,f,sz(2)/2; 0,0,1];
    % compute collineations
    H1 = K * eul([0;p_out(1:2)]) /K;
    H2 = K * eul(p_out(3:5)) /K;
    
end

function [res,J] = costRectif(a_in, sz, m1, m2)
    % costRectif  Compute rectification cost
    
    K3 = commutation(3,3); S = skew([1 0 0]);
    
    n = length(a_in);     % number of variable parameters
    a = [0;0;0;0;0;0];    % default values
    a(1:n) = a_in;        % overwrite default
    
    % intrinsic parameters of the original cameras
    f = 3^a(6)*(sz(1)+sz(2));
    iK = inv([f,0,sz(1)/2; 0,f,sz(2)/2; 0,0,1]);
    % and its Jacobian
    df = (1/f)*log(3);
    DiK = [-df;0;0;0;-df;0;df*sz(1)/2;df*sz(2)/2;0];
    
    % rotation angle X-left is always 0
    [R1,DR1] = eul([0;a(1:2)]);
    [R2,DR2] = eul(a(3:5));
    
    % fundamental matrix btw original points
    F = iK'*R2'*S*R1*iK;
    
    % Jacobian of this parametrization of F
    Jp = [kron(iK', iK'*R2'*S)*DR1(:,2:3),...
        -kron(iK'*R1'*S, iK')*K3*DR2,...
        (-K3*kron(eye(3),iK'*R1'*S*R2) + kron(eye(3),iK'*R2'*S*R1))*DiK];
    Jp(:,n+1:end) = []; % remove fixed parameters
    
    % Sampson residual
    [res, Js] = sampson_fund(F,m1,m2);
    J = Js * Jp;
end


%costRectif  Compute rectification cost
%   a_in is a vector containing the independent variables of the cost
%      function, i.e, five rotation angles  (Y-left, Z-left, X-right,
%      Y-right, Z-right) and (optionally) the focal length. 
%   sz contains size of the image
%   m1, m2 contain the corresponding image points.

% Author: A. Fusiello, 2006, 2018
