function  [G,D] = gpa(M,w)
    %GPA Generalized Procrustes analysis
    % with similarity (7 parameters)
    
    m = length(M);
    
    % constants
    MaxIterations = 50;
    FunctionTol = 1e-3;
    StepTol = 1e-16;
    
    G = mat2cell(repmat(eye(4),1,m),4,4*ones(1,m));
    res_n = cell(1,m); res = Inf;
    for iter = 1: MaxIterations
        
        % compute centroid
        YW = cellfun(@(X,u) bsxfun(@times, u',X), M,w, 'uni',0);
        D = sum(cat(3,YW{:}),3);
        u = 1./sum(cat(3,w{:}),3);  u(isinf(u))=0; % fix divison by 0
        D = bsxfun(@times, u', D);  % same as D = D*diag(u) but faster
        
        % align each M{i} onto centroid
        for i = 1:m
            [R,t,s] = opa(D,M{i},w{i});
            G{i} = [R, t; 0 0 0 1/s] * G{i};
            M{i} = s*(R*M{i} + t);
            res_n{i} = sum(sum((bsxfun(@times,w{i}',(M{i}-D))).^2));
        end
        
        prevres=res;
        res = sum([res_n{:}])/m; % norm(res)
        if norm(res-prevres)<StepTol || norm(res)<FunctionTol
            break; end
        
    end
    
    