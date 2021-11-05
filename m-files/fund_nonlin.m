function  F_out = fund_nonlin(F0, m2, m1)
    %FUND_NONLIN Non-linear refinement of fundamental matrix
 
    % pre-conditioning
    [T1,m1] = precond(m1); [T2,m2] = precond(m2);
    F0 = T2'\(F0/T1);   
    
    % the parametrization is local around U0*V0', 
    % so we have to carry them over
    [f0,U0,V0] = F2par(F0);
    f_out = lsq_nonlin(@(x)fobj(x,m1,m2,U0,V0),f0);
    F_out = par2F(f_out,U0,V0);
    
    %post-conditioning
    F_out = T2' * F_out * T1;
    
end

function [r,J] =  fobj(f,m1,m2,U0,V0)
    % build F from parameters and return residual and its jacobian
    [F,Jp] = par2F(f,U0,V0);
    [r,Js] = sampson_fund(F,m1,m2);
    J = Js * Jp;
end
