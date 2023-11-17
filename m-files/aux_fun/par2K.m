

function [K,J] = par2K( p_in , sz)
%PAR2K Return a K matrix and its jacobian from parameters

n = length(p_in); % number of input parameters
std_width = 2000; % standard image width

if nargin < 2
    sz = [std_width std_width*3/4]; % default size
end
p =  [-1.2, sz(1)/2, sz(2)/2, 1, pi/2]; % default values

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

% Parameters are:  [fov, u0, v0, aspect_ratio, skew]
% The focal lenght is function of hfov and width
% The default hfov is 1.2 rad (68 deg)
%
% sz (optional) is the image size
% used to guess the principal point (u0,v0) when n = 1 (only fov is given)

% K = par2K([fov, u0, v0, ar, skew]) return K with given parameters
% K = par2K([fov, u0, v0, ar]) return K with given parameters and skew=0
% K = par2K([fov, u0, v0]) return K with given parameters and ar=1, skew=0
% K = par2K([fov],sz) return K with given fov and u0=sz(1)/2, v0=sz(2)/2, ar=1, skew=0
% K = par2K([]) returns the K of my mobile phone ;-)


% from 5 to 3 params this is consistent:  par2K(K2par(K,3))
% with 1 parameter only (hfov) the size must be given otherwise the
% principal point is the default

