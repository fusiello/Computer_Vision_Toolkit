function C = imhs(I,s)
%IMHS Harris-Stephens corner strength
% s is the integration scale

    S = fspecial('sobel');
    G = fspecial('gaussian',2*ceil(2*s)+1, s);
     
    % directional derivatives
    Iu = filter2(S, I, 'same');
    Iv = filter2(S',I, 'same');
    
    % convolve with Gaussian
    Iuv = filter2(G, Iu.*Iv,'same');
    Ivv = filter2(G, Iv.^2, 'same');
    Iuu = filter2(G, Iu.^2, 'same');
    
    % trace and determinant
    tr = Iuu + Ivv;
    dt = Iuu.*Ivv - Iuv.^2;
    
    C = dt - 0.04 *tr.^2; % H-S version
    % C = dt./(1+tr);     % Noble version
end
    
    
