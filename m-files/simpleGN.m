function x = simpleGN(res,x0)
    %SIMPLEGN  Gauss-Newton method for non-linear LS
    
    % constants
    MaxIterations = 100;
    FunctionTol = 1e-10;
    StepTol = 1e-6;

    % initialization
    x = x0(:); dx = 0;
    y = Inf;
    
    for iter = 1: MaxIterations
        %% Step 1: Save previous values and take a (tentative) move
        x_prev = x; y_prev = y;
        x = x + dx;
        [y,J]   = res(x);
        fprintf('\tsimpleGN: residual norm %0.5g \n', norm(y));
        
        %% Step 2: Test the move: does it increase error?
        if  norm(y) < norm(y_prev)
            if norm(y-y_prev)<StepTol || norm(y)<FunctionTol
                break; % convergence
            end
        else 
            % restore old values
            x = x_prev; y = y_prev;
            break % GN does not accept error incraese
        end
        
        %% Step 3: Compute a tentative update
        dx = -(J'*J) \ J' * y; % normal equations
        % use pinv for regularizing rank-deficient jacobians
        % dx = -pinv(J) * y;
    end
end





