function [S,J] = skew( u )
% SKEW  Returns the 3x3 skew-symmetric matrix S s.t. skew(u,v) = S*v
% with its jacobian
% Generalized to arbitrary dimensions > 3

u = u(:);
n = length(u);


if n >=3
    
    S = [  0   -u(3)  u(2)
        u(3)    0   -u(1)
        -u(2)  u(1)   0   ];
    
    for i=4:n
        S(end,end+1)=0; % add zero column on the right
        S = [S; [-u(i)*eye(i-1), u(1:i-1)] ]; % add ablock below
    end
    
else
    S = [-u(2) u(1)];
end


if n ==3
    % jacobian is implemented only for n=3
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
    J = [];
    
end
