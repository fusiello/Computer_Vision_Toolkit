function h = H2par(H)
%H2PAR Homography matrix parametrization

H = H./norm(H(3,:));
h = vec(H(1:2,:)');
h(7:8) = cart2ang(H(3,:));
end

