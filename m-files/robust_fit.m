function  [model, inliers] = robust_fit(x, modelfit, modeldist, p, varargin)
    %ROBUST_FIT robust regression
    % Input:
    %      - x: data points
    %      - modelfit: a function that fits the model on the data (with weights)
    %      - modeldist: a function that evaluate the data-model distance
    %      - p: the size of the minimum sample set
    % Options are:
    %      - MSAC: RANSAC improvement (provide inlier threshold)
    %      - LMS:  Least Median of Squares
    %      - IRLS: Iteratively reweighted least squares
    
    opts = parse_options(varargin);
    
    switch opts.method
        case 0
            [model,inliers] = simpleMSAC(x, modelfit, modeldist, p, opts.inlier_th);
        case 1
            [model,inliers] = simpleLMS(x, modelfit, modeldist, p);
        case 2
            [model,inliers] = simpleIRLS(x, modelfit, modeldist);
    end
end

function [bestmodel, bestinliers] = simpleMSAC(x, modelfit, modeldist, p, t)
    %SIMPLEMSAC - Robust fit with the MSAC algorithm
    
    n = size(x,2);
    alpha = 0.99; % Desired probability of success
    f = 0.5 ;     % Pessimistic estimate of inliers fraction
    
    MaxIterations = 1000; % Max number of iterations
    MinIterations = 100;  % Min number of iterations
    mincost =  Inf;
    
    i = 0;
    while  i < max(ceil(log(1-alpha)/log(1-f^p)), MinIterations)
        
        % Generate p random indicies in the range 1..n
        mss = randsample(n, p);
        
        % Fit model to this minimal sample set.
        model = modelfit(x(:,mss),[]);
        
        % Evaluate distances between points and model
        sqres = modeldist(model, x).^2;
        inliers = sqres < (t^2);
        
        % Compute MSAc score
        cost = (sum(sqres(inliers)) + (n -sum(inliers)) * t^2);
        
        if cost < mincost
            mincost = cost;
            bestinliers = inliers; bestmodel = model;
            % Update the estimate of inliers fraction
            f = sum(bestinliers)/n;
        end
        i = i + 1;
        if i > MaxIterations, break, end
    end
    
end


function [bestmodel, bestinliers] = simpleLMS(x, modelfit, modeldist, p)
    %SIMPLELMS - Robust fit with the LMS algorithm
    
    n = size(x,2);
    alpha = 0.99; % Desired probability of success
    f = 0.5 ; % Pessimistic estimate of inliers fraction
    
    MaxIterations  = max( ceil(log(1-alpha)/log(1-f^p)), 100);
    mincost =  Inf;
    
    for  i = 1:MaxIterations
        
        % Generate s random indicies in the range 1..npts
        mss = randsample(n, p);
        % Fit model to this minimal sample set.
        model = modelfit(x(:,mss),[]);
        % Evaluate distances between points and model
        sqres = modeldist(model, x).^2;
        % Compute LMS score
        cost = median(sqres);
        
        scale = 1.4826*sqrt(cost)*(1+5/(length(sqres)-p));
        inliers = sqres < (2.5*scale)^2 ;
        
        if cost < mincost
            mincost = cost;
            bestmodel = model;  bestinliers = inliers;
        end
    end
end

function [model, inliers] = simpleIRLS(x, modelfit, modeldist)
    %SIMPLEIRLS - Robust fit with the IRLS algorithm
    
    MaxIterations = 100;
    FunctionTol = 1e-6;
    StepTol = 1e-6;
    
    % start with OLS fit
    weights =  ones(length(x),1);
    model = modelfit(x, weights);
    obj = Inf;
    
    for  i = 1:MaxIterations
        
        % Evaluate distances between points and model
        res = modeldist(model, x);
        % Estimate scale
        scale = robstd(res);
        % Compute weights
        weights =  weightfun(res,scale);
        % Fit model with weights
        model = modelfit(x, weights);
        
        prevobj= obj;
        obj = sum(weights.*res.^2);
        
        if norm(obj-prevobj)<StepTol || obj<FunctionTol
            break; end
    end
    inliers = abs(res) < (2.5*scale);
    
end


function opts = parse_options(options)
    % parse options and produce an opts struct
    
    % default parameters
    opts.method     = 0;
    opts.inlier_th  = 1.0; % reasonable values btw .5 and 2
    opts.verbose    = false;
    
    % parse optional parameters, if any
    i = 1;
    while i <=  length(options)
        switch options{i}
            case 'MSAC'
                opts.method =  0;
                opts.inlier_th = options{i+1};
                i = i + 1;
            case 'LMS'
                opts.method   = 1;
            case 'IRLS'
                opts.method   = 2;
            case 'Verbose'
                opts.verbose = true;
            otherwise
                error('Robust_fit: unrecognised option');
        end
        i = i + 1;
    end
end














