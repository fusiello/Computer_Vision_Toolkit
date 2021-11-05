function [P_out,M_out,kappa_out] = bundleadj(P0,M0,m,varargin)
    %BUNDLEADJ  Bundle adjustment (Euclidean)
    %
    % Options:
    %  - FixedInterior: Only camera external orientation parameters
    %        and 3D points are adjusted.
    %  - AdjustCommonInterior: Adjust interior parameters, common to
    %       all images. All images have been taken with the same
    %       camera under constant settings.
    %  - AdjustSeparateInterior: Adjust interior parameters, individual
    %       for each image. Each image is captured with a different
    %       camera and/or with varying zoom setting.
    %  - InteriorParameters: When adjusting, specify how many
    %       parameters to use
    %  - FixedPoints: Some (or all) of the 3D points can be kept fixed,
    %       i.e. they can be designated as ground control points.
    % - DistortionCoefficients: use radial distortion either as fixed
    %       values or as initial values to be adjusted (follows the
    %       interior parameters)
    
    assert(size(M0,1) ==3);
    assert(size(m{1},1) ==2);
    
    opts = parse_options(varargin);
    clear varargin;
    
    % no visibility means full visibiity
    if isempty(opts.vis)
        opts.vis  = true(length(M0),length(P0));
    end
    
    % cell of empty arrays is different from an empty cell
    if isempty(opts.kappa)
        opts.kappa  = cell(1,length(P0));
    end
    n_dist_p = length(opts.kappa{1}); % # of distortion parameters 
    
    kappa_out   = cell(1,length(P0));
    
    if opts.verbose
        opts %#ok<NOPRT>
    end
    
    % centralize points (improve conditioning)
    C = mean(M0,2);   M0 = M0 - C;
    Tc = [eye(3), -C];  Tc(4,4) = 1;
    
    n_imm = length(m);
    switch opts.mode
        case 1 % Fixed interior parameters
            n_cam_p = 6; % # of parameters per camera
            index.p = n_cam_p*n_imm;
            index.k = @(i) [];
            index.d = @(i) [];
        case 2  % Adjust common interior parameters (+ radial)
            n_cam_p = 6; % # of parameters per camera
            index.p = n_cam_p*n_imm + opts.n_int_p + n_dist_p;
            index.k = @(i) index.p - (opts.n_int_p + n_dist_p) + 1 : index.p - n_dist_p;
            index.d = @(i) index.p - n_dist_p + 1 :  index.p;
        case 3  % Adjust individual interior parameters (+ radial)
            n_cam_p = 6+opts.n_int_p+n_dist_p; % # of parameters per camera
            index.p = n_cam_p*n_imm;
            index.k = @(i) n_cam_p*(i-1)+7:n_cam_p*(i-1)+7+opts.n_int_p-1;
            index.d = @(i) n_cam_p*(i-1)+7+opts.n_int_p:n_cam_p*i;
    end
    % index.p is end of orientation params
    % index.k are the indices where K is stored 
    % index.d are the indices where distortion coefficients are stored 
    index.r = @(i) n_cam_p*(i-1)+1:n_cam_p*(i-1)+3; % indices of R
    index.t = @(i) n_cam_p*(i-1)+4:n_cam_p*(i-1)+6; % indices of t
    
    % extract initial orientation parameters from the PPMs
    g = zeros(index.p,1);
    K_fix = cell(1,n_imm);
    for i = 1:n_imm
        P0{i} = P0{i} / Tc; % pre-conditioning
        [K,R,t] = krt(P0{i});
        g(index.r(i)) = ieul(R);   g(index.t(i)) = t;
        switch opts.mode
            case 1 % fixed interior parameters.
                K_fix{i} = K;
            case 2  % average common interior parameters (+ radial)
                g(index.k(i)) = g(index.k(i)) + ...
                    (K2par(K,opts.n_int_p)-g(index.k(i)))/i;
                g(index.d(i)) = g(index.d(i)) + ...
                    (opts.kappa{i}-g(index.d(i)))/i;
            case 3  % adjust individual interior parameters (+ radial)
                g(index.k(i)) = K2par(K,opts.n_int_p);
                g(index.d(i)) = opts.kappa{i};     
        end
    end
    
    % partition M in tie-points and control (fixed) points
    M_fix = M0(:,end-opts.gcp+1:end);
    M0 = M0(:,1:end-opts.gcp);
    
    x_out = lsq_nonlin(@(x)fobj(x,m,opts.vis,K_fix,opts.kappa,M_fix,index),[g;M0(:)],'Reduced','Verbose');
    
    % reshape the 3D points (and de-centralize)
    M_out = htx(inv(Tc), [reshape(x_out(index.p+1:end),3,[]), M_fix] );
    
    g_out = x_out(1:index.p);
    % rebuild the PPMs (and de-centralize)
    P_out=cell(1,n_imm);
    for i = 1:n_imm
        switch opts.mode
            case 1 % fixed interior parameters
                P_out{i} = K_fix{i} * [eul(g_out(index.r(i))) g_out(index.t(i))]*Tc;
            case {2,3}  % adjusted interior parameters
                P_out{i} = par2K(g_out(index.k(i))) * [eul(g_out(index.r(i))) g_out(index.t(i))]*Tc;
                kappa_out{i} = g_out(index.d(i));
        end
    end
end

function [r,J,p]  = fobj(x,m,vis,K_fix,kappa_fix,M_fix,index)
    % output: p is the column used to partition J
    
    p = index.p; % lenght of the orientation parameters
    g = x(1:p);
    gcp = size(M_fix,2); % # of fixed points
    M = [reshape(x(p+1:end),3,[]),M_fix];
    
    r = zeros(2*nnz(vis),1);
    J = zeros(2*nnz(vis), p + 3*(size(vis,1)-gcp));
    
    rows = 1:2;
    for i = 1:size(vis,2) % n_imm
        [R,JR] = eul(g(index.r(i)));  t = g(index.t(i));
        [K,JK] = par2K(g(index.k(i)));
        kappa  = g(index.d(i));  F = eye(2);
        
        if isempty(index.k(i))
            K = K_fix{i};
            kappa = kappa_fix{i}; F = double.empty(0,2);
        end
        
        for j = 1: size(vis,1) % n_pts
            if vis(j,i)
                [res,JA,JB,JC,JD] = reproj_res(K, kappa, [R t], M(:,j),m{i}(:,j));
                % same as: J(rows, 6*i-5:6*i) = JA*blkdiag(JR,eye(3));
                J(rows, index.r(i)) = JA(:,1:9) *JR;
                J(rows, index.t(i)) = JA(:,10:end);
                J(rows, index.k(i)) = JC*JK;
                J(rows, index.d(i)) = F*JD;
                
                if  j <= size(vis,1) - gcp % fixed points
                    J(rows, p+3*j-2:p+3*j) = JB;
                end
                r(rows) = res;
                rows = rows+2;
            end
        end
    end
end


function opts = parse_options(options)
    % parse options for BA and produce an opts struct
    
    % initialize default parameters
    opts.mode       = 1;  % fixed interior parameters (do no adjust)
    opts.vis        = []; % visibility (default is full visibility)
    opts.verbose    = false;
    opts.n_int_p    = 3;  % number of interior parameters to use
    opts.kappa      = []; % distortion parameters 
    opts.gcp        = 0;  % number of control points (in the end of M)
    
    % parse optional parameters
    i = 1;
    while i <=  length(options)
        switch options{i}
            case 'FixedInterior'
                opts.mode = 1; % fixed interior parameters (do not adjust)
            case 'AdjustCommonInterior'
                opts.mode = 2; % adjust common interior parameters
            case 'AdjustSeparateInterior'
                opts.mode = 3; % adjust individual interior parameters
            case 'InteriorParameters'
                opts.n_int_p = options{i+1};
                i = i + 1;
            case 'DistortionCoefficients'
                opts.kappa = options{i+1};
                i = i + 1;
            case 'FixedPoints' % use ground control points
                fprintf('\tBA: FixedPoints given\n')
                opts.gcp = options{i+1};% how many gcp
                i = i + 1;
            case 'Visibility'
                fprintf('\tBA: Visibility given\n')
                opts.vis = logical(options{i+1});
                i = i + 1;
            case 'Verbose'
                opts.verbose = true;
            otherwise
                error('\tBA: unrecognised option\n');
        end
        i = i + 1;
    end
    
   % if opts.verbose
        switch opts.mode
            case 1 % Fixed interior parameters
                fprintf('\tBA: FixedInterior\n')
            case 2  % Adjust common interior parameters
                fprintf('\tBA: AdjustCommonInterior\n')
            case 3  % Adjust individual interior parameters
                fprintf('\tBA: AdjustSeparateInterior\n')
        end
        if ~isempty(opts.kappa)
         fprintf('\tBA: Using radial distortion\n')
        end
    %end
    
end

%
%     [~,J] = fobj([g;M(:)],q,ones(size(vis)));
%     figure, spy(J); title('Jacobiana (primaria)'); xlabel('');
%     saveas(gcf, 'jacobiana1', 'epsc')
%     [~,J] = fobj([g;M(:)],q,vis);
%     figure, spy(J); title('Jacobiana (secondaria)'); xlabel('');
%     saveas(gcf, 'jacobiana2', 'epsc')
%     figure, spy(J'*J); title('Hessiana'); xlabel('');
%     saveas(gcf, 'hessiana', 'epsc')
%

