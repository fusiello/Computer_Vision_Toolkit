function [minGSD,maxGSD]  = gsd(K, H, R)
% GSD Ground Sampling Distance (pixel size on the ground)
% K are  the intrinsic parameters
% H is the distance to the ground plane (Z=0), in meters
% R is the rotation matrix of the camera
% minGSD,maxGSD are the min/max GSD, the size in meters of 1 pixel
% projected onto the ground. If R = eye(3) the GSD is constant over the
% image and is equal to H/f, otherwise it varies from point to point.
% This function returns the min and max values over the image.


if nargin < 3
    R = eye(3);
    %same as:
    % xu  = abs(H/K(1,1));
    % xv  = abs(H/K(2,2));
end

a = ieul(R); omega = a(1); phi   = a(2); kappa = a(3);
fu = K(1,1); fv = K(2,2);

% the scale depends on the point
i = 0;
for u =  [-K(1,3), 0 , K(1,3)]
    for v = [- K(2,3), 0, K(2,3)]
        i = i+1;
        xu =  (abs(H)*abs(fu)*abs(fv)*((fv*cos(kappa)*cos(omega) - v*cos(phi)*sin(omega) + fv*sin(kappa)*sin(omega)*sin(phi))^2 + (v*sin(phi) + fv*cos(phi)*sin(kappa))^2)^(1/2))/(fu*fv*cos(omega)*cos(phi) - fu*v*cos(kappa)*sin(omega) + fv*u*sin(kappa)*sin(omega) + fv*u*cos(kappa)*cos(omega)*sin(phi) + fu*v*cos(omega)*sin(kappa)*sin(phi))^2;
        xv =(abs(H)*abs(fu)*abs(fv)*((fu*cos(omega)*sin(kappa) + u*cos(phi)*sin(omega) - fu*cos(kappa)*sin(omega)*sin(phi))^2 + (u*sin(phi) + fu*cos(kappa)*cos(phi))^2)^(1/2))/(fu*fv*cos(omega)*cos(phi) - fu*v*cos(kappa)*sin(omega) + fv*u*sin(kappa)*sin(omega) + fv*u*cos(kappa)*cos(omega)*sin(phi) + fu*v*cos(omega)*sin(kappa)*sin(phi))^2;
        t(i) = norm([xu, xv])/sqrt(2);
    end
end
minGSD = min(t);
maxGSD = max(t);

end



