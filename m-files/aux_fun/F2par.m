function [f,U0,V0] = F2par(F)
    %F2PAR Fundamental matrix parametrization
    
    % the parametrization is local around U0*V0'
    
    [U0,D,V0] = svd(F);
    f(1:6) = [0 0 0 0 0 0];
    f(7) = cart2ang([D(1,1),D(2,2)]);
end

    
    
