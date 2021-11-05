function res  = reproj_res_batch(P,M,m,varargin)
    %REPROJ_RES_BATCH Reprojection residual for all points all cameras
    
    opts = parse_options(varargin);
    clear varargin;
    
    % if no visibility, full visibility
    if isempty(opts.vis)
        opts.vis  = true(length(M),length(P));
    end
    
    if isempty(opts.kappa)
        opts.kappa = num2cell(zeros(1,length(P)),1);
    end
    
    res = [];
    for i = 1:length(P) % n_imm
        
        [K,~,~] = krt(P{i});
        resi =  rdx(opts.kappa{i},htx(P{i}, M), K)  - m{i};
  
        res = [res;resi(repelem(opts.vis(:,i),2)) ];
    end
end

% There are two residual (u and v) per point
% This is the same residual miniized by BA

function opts = parse_options(options)
    % parse options for BA and produce an opts struct
    
    % initialize default parameters
    
    opts.vis        = []; % visibility (default is full visibility)
    opts.kappa      = []; %distortion parameters
    opts.verbose    = false;
    
    % parse optional parameters
    i = 1;
    while i <=  length(options)
        switch options{i}
            case 'DistortionCoefficients'
                opts.kappa = options{i+1};
                i = i + 1;
            case 'Visibility'
                disp('BA: Visibility given')
                opts.vis = logical(options{i+1});
                i = i + 1;
            case 'Verbose'
                opts.verbose = true;
            otherwise
                error('BA: unrecognised option');
        end
        i = i + 1;
    end
end

