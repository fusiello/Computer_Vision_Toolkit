function P = resect_lin(m,M)
%RESECT_LIN  Camera matrix (resection) 
% use DLT algorithm
    
    % pre-conditioning
    [T1,m] = precond(m); 
    [T2,M] = precond(M); % !!! only if M is "shallow" !!!
    
    P = dlt(m, M);
    
    % apply the inverse scaling
    P = T1 \ P * T2;
end













