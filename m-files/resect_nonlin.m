function [P_out] = resect_nonlin(P,m,M)
    %RESECT_NONLIN  Non-linear refinement of resection
    % see also BUNDLEADJ
    
    [K0,R0,t0] = krt(P);
    M = htx ([R0,t0; 0 0 0 1],M); % transform M to work near origin
    
    v_out = lsq_nonlin(@(x)fobj(x,M,m),[0;0;0;0;0;0; K2par(K0)]);
    
    R_out = eul(v_out(1:3));
    P_out = par2K(v_out(7:end)) * [R_out*R0, R_out*t0 + v_out(4:6)];
end


function [res,J]  = fobj(v,M,m)
    
    [K,DK] = par2K(v(7:end)); 
    [R,DR] = eul(v(1:3)); t = v(4:6);
    
    J = []; res = [];
    for i = 1:size(M,2)
        [r,JA, ~, JC] = reproj_res(K,[],[R,t] ,M(:,i),m(:,i));
        J = [J; JA * blkdiag(DR,eye(3)), JC * DK];
        res = [res; r];
    end
end

