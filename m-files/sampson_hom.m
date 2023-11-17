function [rs, Js] = sampson_hom(H,m1,m2)
%SAMPSON_HOM Sampson residual for H (and its Jacobian)

m1=ensure_homogeneous(m1);
m2=ensure_homogeneous(m2);

Q = [1 0; 0 1; 0 0]; [~,DS] = skew([1 1 1]);
K24 = commutation(2,4);
rs = [];  Js = [];
for i = 1:size(m1,2)
    ra = Q'*(skew(m2(:,i))*H*m1(:,i)); % algebraic residual
    
    V = Q'*[skew(m2(:,i))*H*Q, -skew(H*m1(:,i))*Q];
    pV = pinv(V);  % pV = V'*inv(V*V');
    rs = [rs; - pV * ra]; % Sampson residual
    
    if nargout >    1  % Jacobian required
        % derivative of pinv(V)
        DY = - kron(pV',pV) + ( kron((eye(2)-V*pV), pV*pV') + ...
            kron(pV'*pV, (eye(4)-pV*V)) )* K24;
        % derivative of ra
        De = Q'*kron(m1(:,i)',skew(m2(:,i)));
        % derivative of V
        DV = kron(eye(4),Q') * [ kron(Q',skew(m2(:,i)))
            - kron(Q',eye(3)) * DS * kron(m1(:,i)', eye(3)) ];
        
        Js = [Js; - kron(ra',eye(4)) * DY*DV - pV*De];
    end
end
end



% Returns the sampson residual, a 4-vectore whose norm is an
% approximation of the (unsquared) distance from the 4d
% joint-space point [m1;m2] to the H manifold.
%
% Note: there are 4 residuals for each point






