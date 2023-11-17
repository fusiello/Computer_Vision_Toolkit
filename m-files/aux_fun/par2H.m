function [H,J] = par2H(h)
%PAR2H Return a homography matrix and its jacobian from parameters

[v,Jv] = ang2cart(h(7:8));
H = reshape([vec(h(1:6)); v],3,[]).';
J = commutation(3,3)* [eye(6), zeros(6,2); zeros(3,6), Jv];
end