function  H_out = hom_nonlin(H0, m2, m1)
    %HOM_NONLIN Non-linear refinement of H
    
    % pre-conditioning
    [T1,m1] = precond(m1); [T2,m2] = precond(m2);
    H0 = T2*H0/T1; 
    
    h_out = lsq_nonlin(@(x)fobj(x,m1,m2),H2par(H0));
    H_out = par2H(h_out);
    
    % post-conditioning
    H_out = T2\H_out * T1;
end

function [res,J] =  fobj(h,m1,m2)
    % build H from parameters and return residual and its jacobian
    [H,Jp] = par2H(h);
    
    [res, Js] = sampson_hom(H,m1,m2);
    
    J = Js *Jp;
end
