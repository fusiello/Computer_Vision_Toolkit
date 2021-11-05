function [rs, Js] = sampson_hom(H,m1,m2)
    %SAMPSON_HOM Sampson residual for H (and its Jacobian)
    
    m1=ensure_homogeneous(m1);
    m2=ensure_homogeneous(m2);
    
    Q = [1 0; 0 1; 0 0]; [~,DS] = skew([1 1 1]);
    K24 = commutation(2,4); N4 = (eye(4)+commutation(2,2))/2;
    rs = [];  Js = [];
    for i = 1:size(m1,2)
        ra = Q'*(skew(m2(:,i))*H*m1(:,i)); % algebraic residual
        
        V = Q'*[skew(m2(:,i))*H*Q, -skew(H*m1(:,i))*Q]; C = inv(V*V');
        rs = [rs; - V'*C * ra]; % Sampson residual
        
        if nargout >    1  % Jacobian required
            % derivative of V'*(inv(V*V'))
            DY = kron(C,eye(4))*K24 - kron(C,V'*C)*2*N4*kron(V,eye(2));
            % derivative of ra
            De = Q'*kron(m1(:,i)',skew(m2(:,i)));
            % derivative of V
            DV = kron(eye(4),Q') * [ kron(Q',skew(m2(:,i)))
                - kron(Q',eye(3)) * DS * kron(m1(:,i)', eye(3)) ];
            
            Js = [Js; - kron(ra',eye(4)) * DY*DV - V'*C* De];
        end
    end
end



%
% Returns an approximation of the geometric (unsquared) distance
% from the  4d joint-space point [m1;m2] to the H manifold.
%
% Note: there are 4 residuals for each point







