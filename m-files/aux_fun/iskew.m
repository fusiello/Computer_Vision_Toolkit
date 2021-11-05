function  x = iskew(S)
%ISKEW Returns the  vector x s.t. S = skew(x)  
%
%   S is  skew-symmetric matrix  s.t. cross(x,y) = S*y
%   

% Author: A. Fusiello


if (size(S,1) ~=3) ||  (size(S,2) ~=3)
    error('Argument must be a 3x3 matrix');
end


%x=[S(3,2) S(1,3) S(2,1)]';

x= [(S(3,2)-S(2,3))/2 (S(1,3)-S(3,1))/2 (S(2,1)-S(1,2))/2]';
