function [x,J]  = lsq_nonlin(res,x0, varargin)
    %LSQ_NONLIN  Levenberg-Marquardt method for non-linear LS
    %
    % Options:
    %   'GaussNewton'    - solve with Gauss-Newton (no damping)
    %   'Reduced'        - solve reduced normal equations (only for BA)
    
    use_LSQNONLIN = false;
    
    % If the optimization toolbox is not installed, use replecememt
    if ~license('test','optimization_toolbox') 
        use_LSQNONLIN = false;
    end
    
    opts = parse_options(varargin);
    
    % constants
    MaxIterations = 200;
    FunctionTol = 1e-10;
    StepTol = 1e-6;
    
    lambda = opts.lambda0;
    % lambda = 0 is Gauss-Newton
    % lambda > 0 is Levenberg-Marquardt (default)
    
    if ~use_LSQNONLIN 
        % replacement for lsqnonlin
        
        % initialization
        x = x0(:); dx = 0;
        y = Inf;
        J = [];
        
        for iter = 1: MaxIterations
            %% Step 1: Save previous values and take a (tentative) move
            x_prev = x; y_prev = y;  J_prev = J;
            
            x = x + dx;
            if opts.reduced
                [y,J,p] = res(x);
                % p is the column used to partition J
                J = sparse(J);
            else
                [y,J]   = res(x);
            end
            % here we have new value for x,y,J
            
            if opts.verbose && iter == 1
                fprintf('\tlsq_nonlin: initial residual norm %0.5g \n', norm(y));
            end
            
            %% Step 2: Test the move: does it increase error?
            
            if  norm(y) < norm(y_prev)
                % if the new error is less than the old error, then
                % accept the new values of the parameters, diminish
                % the value of lambda by a factor of 10
                lambda = lambda/10;
                
                % termination test (only if error decrease)
                if norm(y-y_prev)<StepTol || norm(y)<FunctionTol
                    break;
                end
            elseif iter>1
                % at first iteration prev values are meaningless
               
                % if the new error is greater than the old error, then
                % keep the old parameter values, increase the value of
                % lambda by a factor of 10
                lambda = lambda*10;
                x = x_prev;
                y = y_prev;
                J = J_prev;
                
                % GN does not accept error incraese: exit
                if lambda==0, break, end
                
            end
            % The move is effective or the previous state is restored
            
            %% Step 3: Compute a tentative update
            
            H = J'*J; % Hessian
            % LM damping + fix zero entries on the diagonal of H
            D = diag(diag(H)  + (abs(diag(H))<eps) );
            H = H + lambda*D;  % damped Hessian
            
            if opts.reduced % only for BA
                % solve reduced normal equations of BA
                b = -J'*y;
                Hpc = H(p+1:end,1:p);
                iHpp = H(p+1:end,p+1:end) \ speye(size(Hpc,1)); % inv(Hpp)
                Htmp =  Hpc'*iHpp;  % Hcp*iHpp;
                Hr = (H(1:p,1:p)-Htmp*Hpc );
                % jacoby preconditioner 1/diag(Hr)
                P = spdiags(1./diag(Hr),0,size(Hr,1),size(Hr,1));
                
                dxc = (P*Hr) \ (P*(b(1:p)-Htmp*b(p+1:end)));
                dxp = iHpp* (b(p+1:end) - Hpc*dxc) ;
                dx = [dxc;dxp];
            else
                dx = - H \ J' * y;
            end
            % here we have a tentative update dx
        end
        if opts.verbose 
            fprintf('\tlsq_nonlin: final residual norm %0.5g \n', norm(y));
        end
    else % use lsqnonlin from optimization toolbox 
        
        options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
            'SpecifyObjectiveGradient',true,'MaxIterations',MaxIterations,...
            'StepTolerance', StepTol, 'FunctionTolerance', FunctionTol,...
            'CheckGradients',true, 'Display','off', 'FiniteDifferenceType','central');
        
        [x,~,~,~,~,~,J] = lsqnonlin(res,x0(:),[],[],options);
    end
end


function opts = parse_options(options)
    % parse options for LM and produce an opts struct
    
    % default parameters
    opts.reduced   = false;
    opts.lambda0   = 1;
    opts.verbose   = false;
    
    % parse optional parameters, if any
    i = 1;
    while i <=  length(options)
        switch options{i}
            case 'Reduced'
                opts.reduced   = true;
            case 'GaussNewton'
                opts.lambda0 =  0;
            case 'Verbose'
                opts.verbose = true;
            otherwise
                error('lsq_nonlin: unrecognised option');
        end
        i = i + 1;
    end
end



