function K_out = autocal(F,K0)
    %AUTOCAL Autocalibration from fundamental matrices
    
    k0 = K2par(K0,3);  % aspect ratio=1, skew=0
    
    k_out = lsq_nonlin(@(x)fobj(x,F),k0);
    
    K_out = par2K(k_out);

end

function [res,J]  = fobj(a,F)
    
    J =[]; res=[];
    for i=1:size(F,1)
        for j=1:size(F,2)
            if ~isempty(F{i,j})
                [r, D] = jacobianHF(a,F{i,j});
                res = [res; r];
                J = [J; D];
            end
        end
    end
end


function [r,J] = jacobianHF(k,F)
%JACOBIANHF  Jacobian of the Huang-Faugeras residual
       
    K3 = commutation(3,3);    
    [K,JK] = par2K(k);
    
    E = K'*F*K; EE = E'*E;
    c1 = 2*trace(EE*EE);
    c2 = (trace(EE))^2;
    
    r = c1/c2-1; % objective f.
    
    JH = [ 1/c2, -c1/(c2^2)];
    JS = 4* [ EE(:)' * (eye(9)+K3) * kron(eye(3),E')
              trace(EE)*(E(:))'];
    JE = K3*kron(eye(3),K'*F') + kron(eye(3),K'*F);
    J = JH * JS * JE * JK; % jacobian
end
