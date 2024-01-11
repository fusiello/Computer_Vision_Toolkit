function p = K2par(K, n, sz )
%K2PAR K matrix parametrization
% n is the number of parameters required (1 to 5)
% sz is the image size

% sz must be provided for p to be meaningful, 
% however par2K(K2par(K, 5, sz), sz) is consistent with any sz, 
% so we provide a default only for this usage
if nargin < 3 
    sz  = [2000, 1500];
end

if nargin < 2
    n = 5; % default is 5 parameters
end

hfov = 2*atan(sz(1)/(2*K(1,1)));
theta = atan(-K(1,1)/K(1,2));

p = [hfov,  K(1,3), K(2,3),  (K(2,2)/K(1,1))*sin(theta),  theta];
p = p(1:n)'; % clip to the number of reqired parameters
end

% Parameters are:  [hfov, u0, v0, ar, skew]
% The focal lenght is function of hfov and width.





