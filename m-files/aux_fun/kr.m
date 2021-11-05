function X = kr(A,B)
%KR Khatri-Rao product.
%   kr(A,B) returns the Khatri-Rao product of two matrices A and B, of 
%   dimensions I-by-K and J-by-K respectively. The result is an I*J-by-K
%   matrix formed by the matching columnwise Kronecker products, i.e.
%   the k-th column of the Khatri-Rao product is defined as
%   kron(A(:,k),B(:,k)).
%
%   If either X or Y is sparse, only nonzero elements are multiplied
%   in the computation, and the result is sparse.
%
%   See also kron.

%   Authors: Andrea Fusiello, 2016 based on Laurent Sorber's kr and Matlab
%   kron


%A = sparse(A);
%B = sparse(B);

[ma,na] = size(A);
[mb,nb] = size(B);

if (na ~=nb)
    error('kr:ColumnMismatch', ...
        'Input matrices must have the same number of columns.');
end


if ~issparse(A) && ~issparse(B)
    
    % Both inputs full, result is full.
   
    A = reshape(A,[1 ma na]);
    B = reshape(B,[mb 1 na]);
    X = reshape(bsxfun(@times,B,A),[ma*mb  na]);

else
    
    % At least one input is sparse, result is sparse.

    val=[]; I=[]; J=[];
    
    for k = 1:size(A,2)
        
        %X(:,k) = kron(A(:,k),B(:,k));
        [ia,~,sa] = find(A(:,k));
        [ib,~,sb] = find(B(:,k));
        
        ia = ia(:);  sa = sa(:);
        ib = ib(:);  sb = sb(:);
        
        ik = bsxfun(@plus, mb*(ia-1).', ib);
        jk = ones(size(ik));
        
        s = bsxfun(@times,sb,sa.');
        
        val = [val;s(:)];
        I = [I;ik(:)];
        J = [J;k*jk(:)];

    end
    
    X = sparse(I,J,val,ma*mb,na) ;
    
end
    
