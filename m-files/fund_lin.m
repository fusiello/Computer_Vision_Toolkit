function F = fund_lin(m2,m1,w)
%FUND_LIN  Fundamental matrix with 8-points algorithm

if nargin < 3 || isempty(w)
    w = ones(size(m1,2),1);  % weights
end

% pre-conditioning
[T1,m1] = precond(m1);
[T2,m2] = precond(m2);

F = eight_points(m2, m1, w);

% enforce singularity of  F
[U,D,V] = svd(F);
D(3,3) = 0; F = U*D*V';

% apply the inverse scaling
F = T2' * F * T1;

end

