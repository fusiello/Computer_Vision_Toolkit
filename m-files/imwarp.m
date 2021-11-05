function [Iout, mask]  = imwarp(I,f,bb)
%IMWARP  Image warp 
% Apply the projective transformation specified by the _inverse_ of f 
% The bounding box is specified with [minx; miny; maxx; maxy];
% mask is 1 for valid pixels and 0 otherwise

[x,y] = meshgrid(bb(1):bb(3),bb(2):bb(4));
% pp = htx(inv(H),[x(:),y(:)]');
pp = f([x(:),y(:)]');

xi=reshape(pp(1,:),size(x,1),[]);
yi=reshape(pp(2,:),size(y,1),[]);
Iout=interp2(1:size(I,2),1:size(I,1),double(I),xi,yi,'linear',NaN);

mask = (~isnan(Iout));
Iout(~mask) = median(I(:));  % the background value is the median of the image

cast(Iout,class(I)); % cast I2 to whatever was I

% implements backward mapping, so f is the inverse of the function
% that is being applied to I
