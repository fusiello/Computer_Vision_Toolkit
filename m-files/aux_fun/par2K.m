function [K,J] = par2K( p_in , sz)
%PAR2K Return a K matrix and its jacobian from parameters
% sz is the image size

n = length(p_in); % number of input parameters

% sz must be provided for K to be meaningful, 
% however par2K(K2par(K,5,sz), sz) is consistent with any sz, 
% so we provide a default only for this usage
if nargin < 2
    sz  = [2000, 1500];
end

% default parameters
p =  [-1.2, sz(1)/2, sz(2)/2, 1, pi/2];
% the dafault fov is 1.2 rad = 68 deg

p(1:n) = p_in;   % overwrite default

a =  0.5*sz(1)/tan(0.5*p(1));

K = [ a, -a/tan(p(5)),       p(2)
      0,  a*p(4)/sin(p(5)),  p(3)
      0   0                  1   ];

J = zeros(9,5);
J(1)= (0.5*sz(1))/(cos(p(1)) - 1);
J(4)=(0.25*sz(1)*(tan(0.5*p(1))^2 + 1))/(tan(0.5*p(1))^2*tan(p(5)));
J(5)= -(0.25*p(4)*sz(1))/(sin(0.5*p(1))^2*sin(p(5)));
J(16)=1;
J(26)=1;
J(32)=(0.5*sz(1))/(tan(0.5*p(1))*sin(p(5)));
J(40)=(0.5*sz(1)*(tan(p(5))^2 + 1))/(tan(0.5*p(1))*tan(p(5))^2);
J(41)=-(0.5*p(4)*sz(1)*cos(0.5*p(1))*cos(p(5)))/(sin(0.5*p(1))*sin(p(5))^2);

J(:,n+1:end) = []; % cut the fixed parameters

end

% Parameters are:  [hfov, u0, v0, aspect_ratio, skew]
% sz is used to compute hfov and to guess (u0,v0) if needed
% K = par2K([hfov, u0, v0, ar, skew], sz) return K with given parameters
% K = par2K([hfov, u0, v0, ar], sz) return K with given parameters and skew=0
% K = par2K([hfov, u0, v0], sz ) return K with given parameters and ar=1, skew=0
% K = par2K([hfov],sz) return K with given hfov and u0=sz(1)/2, v0=sz(2)/2, ar=1, skew=0
% K = par2K([], sz) returns the K with hfov=1.2, u0=sz(1)/2, v0=sz(2)/2, ar=1, skew=0


