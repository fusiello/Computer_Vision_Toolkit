function E = eight_points(m2, m1, w)
%EIGHT_POINTS 8-points algorithm for F/E

    if nargin <3, w = ones(size(m1,2),1); end % weights

    m1 = ensure_homogeneous(m1);  
    m2 = ensure_homogeneous(m2);
    
    L = [];
    for i = 1: size(m1,2)
        L = [L; w(i) * kron( m1(:,i)' , m2(:,i)') ];
    end
    
    [~,~,V] = svd(L);
    E = reshape(V(:,end),3,3);
    
end
