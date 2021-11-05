function [R,t,s] = opa(D,M,w)
    %OPA Orthogonal Procrustes analysis
    % Solves the othogonal procrustes problem D=s*(R*M+t)
    
    if nargin==3  % Weights  provided
        w=w(:); % make sure it is a column
        % same as  D = D*W;  M = M*W;  but faster
        D = bsxfun(@times, w', D);
        M = bsxfun(@times, w', M);
        
    else  % otherwise no weights
        w = ones(size(M,2),1);
    end
    
    J = eye(size(M,2)) - (w*w')/(w'*w);  % centering matrix
    
    [U,S,V] = svd(D*J*M'); % solve for rotation
    R = U*diag([1,1,det(V'*U)])*V';
    
    if nargout==2  % no scale
        s = 1;
    else   % solve for scale
        s = trace(S) / trace(M*J*M');
        % trace(D*J*M'*R') == trace(S)
    end
    
    t =(D/s -R*M)*w/(w'*w); % solve for translation
end
