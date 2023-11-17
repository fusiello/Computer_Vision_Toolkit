function p = K2par( K, n, std_width)
%K2PAR K matrix parametrization
% n is the number of parameters required
% std_width is the standard widh used to compute the fov

if nargin < 3
    std_width = 2000;
end

if nargin < 2
    n = 5; % default is 5 parameters
end

hfov = 2*atan(std_width/(2*K(1,1)));
theta = atan(-K(1,1)/K(1,2));

p = [hfov,  K(1,3), K(2,3),  (K(2,2)/K(1,1))*sin(theta),  theta];
p = p(1:n)'; % clip to the number of reqired parameters
end

% Parameters are:  [hfov, u0, v0, ar, skew]
% The focal lenght is function of hfov and width.
% By default width is set to a standard value of 2000.
% This value is not  important as long as the same value is used also in
% par2K, but it is understood that hfov is fictious. If one wants the
% *real* hfov he/she must supply the *real* width.



