function   [R,t]  = relative_nonlin(R0, t0, q2, q1, K2, K1)
    %RELATIVE_NONLIN Non-linear refinement of relative orientation
   
    q1 = K1\ensure_homogeneous(q1); 
    q2 = K2\ensure_homogeneous(q2); 
    % now q1 and q2 are homogeneous NIC
    
    q1 = R0*q1; % tranform q1 to work near origin
    
    e_out = lsq_nonlin(@(x)fobj(x,q1,q2),[0;0;0;cart2ang(t0)]);
    
    R = eul(e_out(1:3))*R0;  t = ang2cart(e_out(4:5));
end

function [res,J] =  fobj(e,q1,q2)
    % build F from parameters and return residual and its jacobian
    [E,Jp] = par2E(e);
    
    [res,Js] = sampson_fund(E,q1,q2);
    
    J = Js *Jp;
end
