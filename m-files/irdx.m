function m_i = irdx(k_in,m_d,K)
    %IRDX Undo radial distortion
    
    n = length(k_in);   % # of parameters
    kappa = zeros(4,1); % default values (up to 4)
    kappa(1:n) = k_in;  % overwrite default
    
    % convert to normalizaed image coordinates
    m_d = htx(inv(K),m_d);
    
    % for each point independently
    for i=1:size(m_d,2)
        m_i(:,i) = lsq_nonlin(@(x)fobj(x,m_d(:,i),kappa),m_d(:,i));
    end
    
    % convert back to pixels
    m_i = htx(K,m_i);
end

function [res,J]  = fobj(m_ideal,m_dist,k)
    
    [m,J] = rdx(k,m_ideal,eye(3)); % in NIC
    res = m - m_dist;
end

% input coordinates are distorted
% k up to 4 coefficients
