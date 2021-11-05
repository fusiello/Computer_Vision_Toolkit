function k = radial(m_dist, m_ideal, K,n)
    %RADIAL Regress radial distortion parameters 
    
    % convert to normalizaed image coordinates
    m_dist = htx(inv(K),m_dist);
    m_ideal = htx(inv(K),m_ideal);
    
    L=[];
    r2 = ( m_ideal(1,:).^2 + m_ideal(2,:).^2 ); % square radius
    for i = 1:size(m_ideal,2)    
        L = [L; kron(ones(1,n), m_ideal(:,i)) .* (r2(i).^(1:n)) ];
    end
    % linear regression
    k = L \ (m_dist(:)-m_ideal(:));
end

 % k is s.t. to warp m_ideal on m_dist
 % n is the the number of distortion  parameters 
