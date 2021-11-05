function A = blocksum(X,n)
    %  BLOCKSUM  sums the square blocks of order n
    
    % The dimensions of the input matrix A must be multiple of sublen
    
    if rem(size(X,1),n)
        error('X must have size multiple of n')
    end
    
    nnodes = size(X,1)/n;
      
    if issparse(X)
        D = kron(speye(nnodes), ones(1,n));
        A = (D*X*D');
    else
        % equivaent to above but faster for full matrices
        A = squeeze(sum(reshape(sum(reshape(X,n,[])),size(X,1)/n,n,[]),2));
        
    end
    
end

