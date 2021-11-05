function [S,J] = skew( u )
    % SKEW  Returns the skew-symmetric matrix S s.t. skew(u,v) = S*v
    % Also return the jacobian
    % Generalized to arbitrary dimensions
    
    u = u(:);
    n = length(u);
    
    if n == 3
        
        S = [   0    -u(3)  u(2)
            u(3)    0   -u(1)
            -u(2)  u(1)   0   ];
        
        J = [0     0     0
             0     0     1
             0    -1     0
             0     0    -1
             0     0     0
             1     0     0
             0     1     0
            -1     0     0
             0     0     0];
    else
        
        S =  [];
        for i=1:n
            S = [S;
                zeros(n-i, i-1), -u(i+1:end), u(i) * eye(n-i)];
        end
        
        J = []; % not implemented
        
    
        % if n == 3
        %     S = [0 0 1; 0 -1 0; 1 0 0] * S;
        % end
        
    end
    
    % if n==3 S is the skew-symmetric matrix [u]_x up to:
    % swapping 1st and 3rd rows and  changing the sign of the 2nd row
    % S=[   0    -x(3)  x(2)
    %      x(3)    0   -x(1)
    %     -x(2)  x(1)   0   ];
