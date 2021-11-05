function M = triang_nonlin(M0, P, m)
    %TRIANG_NONLIN Non-linear refinement of triangulation
    
    M = lsq_nonlin(@(x)fobj(x,P,m),M0,'GaussNewton');
end

function [res,J] = fobj(M,P,m)
    J=[];  res=[];
    for j=1:length(P)
        [K,R,t ] = krt(P{j});
        [r,~,JB] = reproj_res(K, [], [R,t], M, m{j});
        J = [J ; JB];
        res = [res; r];
    end
end





