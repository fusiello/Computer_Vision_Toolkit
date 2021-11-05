function [R,t] = exterior_nonlin(R0,t0,q,M,K)
    %EXTERIOR_NONLIN Non-linear refinement of exterior orientation
    
    q = K\ensure_homogeneous(q);  % Homogeneous NIC
    
    M = htx ([R0,t0; 0 0 0 1],M); % transform M to work near origin
    
    g_out = lsq_nonlin(@(x)fobj(x,M,q(1:2,:)),[0;0;0;0;0;0]);
    
    R_out = eul(g_out(1:3));
    R =  R_out*R0;  t =  R_out*t0 + g_out(4:6);
end

function [res,J]  = fobj(g,M,q)
    
    [R,JR] = eul(g(1:3));  t = g(4:6);
    
    J =[];  res = [];
    for i = 1:size(M,2)
        [r,JA] = reproj_res(eye(3),[],[R,t], M(:,i),q(:,i));
        J = [J ; JA * blkdiag(JR,eye(3))];
        res = [res; r];
    end
end
