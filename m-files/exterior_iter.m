function [R,t] = exterior_iter(m,M,K)
 %EXTERIOR_ITER Exterior orientation with anisotropic Procrustes analysis  
 
    Q = K\ensure_homogeneous(m);   % homogeneous NIC

    MaxIterations = 500;
    FunctionTol = 1e-6;
    StepTol = 1e-6;
    
    obj = +Inf;
    z = ones(1,size(Q,2));
    for  i = 1:MaxIterations
        
        % solve for R,t
        [R,t] = opa(Q*diag(z),M);
        B = R*M + t;

        % update objective and test convergence
        prevobj= obj;
        obj = norm(B - Q*diag(z),'fro');      
        if norm(obj-prevobj)<StepTol || obj<FunctionTol
            break; end
        
        % solve for z
        z  =(kr(eye(size(Q,2)),Q))\B(:);
        %  fix sign of scale
        z = z * sign(z(1));
        z(z<0) = 0;
    end
  
    % Input
    % m : matrix (3xn) image coordinates in image-space [u v 1]
    % M : matrix (3xn) coordinates in object-space [X Y Z]
    % Output
    % R : matrix (3x3) rotation (world-to-camera)
    % t : vector (3x1) translation (world-to-camera)
    
    % Algorithm ref.: Garro, Crosilla, Fusiello, 3DIMPVT 2012
    
