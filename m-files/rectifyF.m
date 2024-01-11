function [H1,H2, K] = rectifyF(m1, m2, sz)
% RECTIFYF Epipolar rectification (uncalibrated)

m1=ensure_homogeneous(m1);
m2=ensure_homogeneous(m2);

% first optimize with fixed f (to avoid jacobian singularity)
p_init = lsq_nonlin(@(x)costRectif(x,sz,m1,m2), [0 0 0 0 0]);

% then add f to the optimization
p_out = lsq_nonlin(@(x)costRectif(x,sz,m1,m2), [p_init(:); -1.2]);

% if extreme focal, fall back to default one
if abs(p_out(6)) <.26 || abs(p_out(6)) > 2.3
    p_out = [p_init(:); -1.2];
end

% recover K
K = par2K(p_out(6), sz);

% compute homographies
H1 = K * eul([0;p_out(1:2)]) /K;
H2 = K * eul(p_out(3:5)) /K;

end

function [res,J] = costRectif(a_in, sz, m1, m2)
% costRectif  Compute rectification cost

K3 = commutation(3,3); S = skew([1 0 0]);

n = length(a_in);     % number of variable parameters
a = [0;0;0;0;0;-1.2]; % default values
a(1:n) = a_in;        % overwrite default

% intrinsic parameters of the original cameras
[K, JK] = par2K(a(6),sz);
iK = inv(K);

% and its Jacobian
DiK = -kron(iK', iK) * JK;

% rotation angle X-left is always 0
[R1,DR1] = eul([0;a(1:2)]);
[R2,DR2] = eul(a(3:5));

% fundamental matrix btw original points
F = iK'*R2'*S*R1*iK;

% Sampson residual
[res, Js] = sampson_fund(F,m1,m2);

% Jacobian of this parametrization of F
Jp = [kron(iK', iK'*R2'*S)*DR1(:,2:3),...
    -kron(iK'*R1'*S, iK')*K3*DR2,...
    (-K3*kron(eye(3),iK'*R1'*S*R2) + kron(eye(3),iK'*R2'*S*R1))*DiK];
Jp(:,n+1:end) = []; % remove fixed parameters

J = Js * Jp;
end


%costRectif  Compute rectification cost
%   a_in is a vector containing the independent variables of the cost
%      function, i.e, five rotation angles  (Y-left, Z-left, X-right,
%      Y-right, Z-right) and (optionally) the focal length.
%   sz contains size of the image
%   m1, m2 contain the corresponding image points.

% Author: A. Fusiello, 2006, 2018

% The focal lenght is parametrized with the horizontal fov
% The default hfov is 1.2 rad (68 deg) 
