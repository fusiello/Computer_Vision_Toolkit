function  [rs, Js] = sampson_fund(F,m1,m2)
    %SAMPSON_FUND Sampson residual for F (and its Jacobian)
    
    m1=ensure_homogeneous(m1);
    m2=ensure_homogeneous(m2);
    
    K3 = commutation(3,3); Q = [1 0; 0 1; 0 0];
    rs = [];  Js = [];
    for i = 1:size(m1,2)
       ra = m2(:,i)'*F*m1(:,i); % algebraic residual
        
        v = [ m2(:,i)'*F*Q , m1(:,i)'*F'*Q ];   c = 1/(v*v');
        rs = [rs; -v'*c * ra];   % Sampson residual
        
        if nargout > 1 % Jacobian required
            %  derivative of J'*(inv(J*J'))
            DY = (c*eye(4) - 2*c^2 *(v'*v));
            % derivative of ra
            De =  kron(m1(:,i)', m2(:,i)');
            % derivative of v
            Dv = [kron(Q', m2(:,i)'); kron(Q', m1(:,i)')*K3];
            
            Js =  [Js;  - ra*DY*Dv - v'*c*De];
        end
    end
end


% Returns the sampson resisual, a 4-vectore whose norm is an
% approximation of the distance from the 4d joint-space point
% [m1;m2] to the F manifold.
%
% Note: there are 4 residuals for each point
