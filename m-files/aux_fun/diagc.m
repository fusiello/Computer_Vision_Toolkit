function B = diagc(A)
%%DIAGC  Block diagonal costruction.
%   columns of A becomes the diagonal blocks of B

C = num2cell(A,1);
B = blkdiag(C{:});

% equivalent to:
% [m,n] = size(A);
% B =kron(eye(n),ones(m,1)) .* kron(ones(n,1), A);

end

