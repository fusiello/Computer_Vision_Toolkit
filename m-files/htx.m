function Y = htx( T, X )
%HTX Apply homogeneous transform 
%   Apply transformation T in homogeneous coordinates
%   T can be any dimension 
%   T 3x4 is a projection
%   T 3x3 is a transformation of the projective plane
%   T 4x4 is a transformation of the projective space

Y = T * [X; ones(1, size(X,2))]; % apply tranform

nr = size(Y,1);

Y = Y ./ repmat(Y(nr,:),nr,1); % perspective division

Y(nr,:) = []; % remove trailing ones

end
