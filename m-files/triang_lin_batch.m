function X = triang_lin_batch(P, x, vis)
    %TRIANG_LIN_BATCH  Triangulation for n points in multiple images
    
    % if no visibility, full visibility
     if nargin < 3
        vis  = true(size(x{1},2),length(P));
    end
    
    X=NaN(3,size(x{1},2));
    for i = 1:size(x{1},2)
        % triangulationfor the i-th point
        L=[];
        % stack equations for every camera
        for j=1:length(P)
            if vis(i,j)
                L = [L; skew([x{j}(:,i);1])*P{j} ]  ;
            end
        end
        [~,~,V] = svd(L);
        X(:,i)  = V(1:3,end)./V(4,end);
    end
end


