function P = dlt( x, X, w)
%DLT Direct Linear Transform
    
   if nargin <3, w = ones(size(X,2),1); end % weights

    X = [X; ones(1, size(X,2))];
    x = [x; ones(1, size(x,2))];
    
    L = [];
    for i = 1: size(X,2)
        L = [L; w(i) * kron( X(:,i)', skew(x(:,i)) )];
    end
    
    [~,~,V] = svd(L);
    P = reshape(V(:,end),size(x,1),[]);
end

    
