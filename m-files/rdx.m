function [m_d,Jm,Jk] = rdx(k_in,m_i,K)
    %RDX Apply radial distortion
    
    n = length(k_in); % # of parameters
    k = zeros(4,1);  % default values (up to 4)
    k(1:n) = k_in;    % overwrite default
    
    % convert to normalizaed image coordinates
    m_i = htx(inv(K),m_i);
    x = m_i(1,1); y = m_i(2,1);
    r2 = ( m_i(1,:).^2 + m_i(2,:).^2 ); % square radius
    
    m_d = m_i.*(k(1)*r2 + k(2)*r2.^2 + k(3)*r2.^3 + k(4)*r2.^4 + 1);
    
    % convert back to pixels
    m_d = htx(K,m_d);
    
    % Jacobian matrices for a **single** point
    Jm = K(1:2,1:2) * ...
        [x*(2*k(1)*x + 4*k(2)*x*r2(1) + 6*k(3)*x*r2(1)^2 + 8*k(4)*x*r2(1)^3)...
         + k(2)*r2(1)^2 + k(3)*r2(1)^3 + k(4)*r2(1)^4 + k(1)*r2(1) + 1,...
         x*(2*k(1)*y + 4*k(2)*y*r2(1) + 6*k(3)*y*r2(1)^2 + 8*k(4)*y*r2(1)^3);
         y*(2*k(1)*x + 4*k(2)*x*r2(1) + 6*k(3)*x*r2(1)^2 + 8*k(4)*x*r2(1)^3),...
         y*(2*k(1)*y + 4*k(2)*y*r2(1) + 6*k(3)*y*r2(1)^2 + 8*k(4)*y*r2(1)^3)...
         + k(2)*r2(1)^2 + k(3)*r2(1)^3 + k(4)*r2(1)^4 + k(1)*r2(1) + 1]...
         * inv(K(1:2,1:2)); %#ok<MINV>
    
    Jk = double.empty(2,0);
    for i = 1:n
        Jk = [Jk, K(1:2,1:2) * [x;y].* r2(1).^i];
    end
    
end

% output distorted coordinates according to the coefficients in k_in
% up to 4 coefficients are considered
% The Jacobian of m_d wrt k_in is also returned (for BA)\

