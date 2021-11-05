function [k,P_est] = calib_radial( x_dist, X, n)
    %CALIB_RADIAL calibrate radial distortion with resection
    %  n is the number od coefficient
    
    % constants
    StepTol = 1e-6;
    MaxIterations = 100;
    
    res = +Inf;
    x_u = x_dist; % x_u are the undistorted coordinates
    k = 0;
    for iter = 1: MaxIterations
        
        % resection with current undistorted
        P_est = resect_lin(x_u, X);
        [K,~,~]  =  krt(P_est);
        
        % project with curent PPM
        x_est = htx(P_est,X);
        
        % compute radial distortion
        k = radial(x_dist, x_est, K, n);
        
        % update undistorted with current estimate of distortion
        x_u = irdx(k,x_dist,K);
        
        % convergence test
        prev_res = res;
        res =  rmse(x_u(:)-x_est(:));
        if abs(res-prev_res)<StepTol
            break;
        end
    end
end

% Warning: this procedure is rather unstable, it can converge to wrong
% result depending on the distribution of points in the image



