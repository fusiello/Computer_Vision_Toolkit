function [K,J] = par2K( p_in )
    %PAR2K Return a K matrix and its jacobian from parameters
    
    n = length(p_in); % number of parameters
    p =  [0,0,1,1,0]; % default values
    p(1:n) = p_in;    % overwrite default
    
    f = p(3)*( p(1)+p(2) );
    
    K = [ f,    p(5),    p(1)
        0,     f*p(4),   p(2)
        0      0          1 ];   
    
    J = [  p(3),  p(3),     p(1) + p(2),            0,            0
           0,     0,            0,                  0,            0
           0,     0,            0,                  0,            0
           0,     0,            0,                  0,            1
        p(3)*p(4), p(3)*p(4), p(4)*(p(1)+p(2)), p(3)*(p(1)+p(2)), 0
           0,     0,            0,                  0,            0
           1,     0,            0,                  0,            0
           0,     1,            0,                  0,            0
           0,     0,            0,                  0,            0];
    
    J(:,n+1:end) = []; % cut the fixed parameters
end




% Parameters are:  [u0, v0, focal, aspect_ratio, skew]
% The focal lenght is parametrized wrt (u0,v0) 
% The default focal length is u0+v0, which gives a diagonal view
% angle of approx. 60 deg
% Plese note that the first two default values are meaningless,
% so DO NOT use this with less than 2 parameters
