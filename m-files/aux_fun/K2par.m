function p = K2par( K, n )
    %K2PAR K matrix parametrization
    % n is the number of parameters required 
    
    if nargin<2
        n = 5; % default is 5 parameters
    end
    
    p = [K(1,3),K(2,3),K(1,1)/(K(1,3)+K(2,3)),K(2,2)/K(1,1),K(1,2)];    
    p = p(1:n)'; % clip to the number of reqired parameters
end


% Parameters are:  [u0, v0, focal, aspect_ratio, skew]
% The focal lenght is parametrized wrt (u0,v0) 
% DO NOT use this with less than 2 parameters 


