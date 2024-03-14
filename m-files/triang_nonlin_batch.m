function M =triang_nonlin_batch(M0, P, m, vis)
%TRIANG_NONLIN_BATCH  Non-linear refinement of triangulation for n points
%in multiple images

% if no visibility, full visibility
if nargin < 4
    vis  = true(size(m{1},2),length(P));
end

 % reshape the points cell matrix
    points  = cat(3,m{:});
    
M = NaN(3,size(m{1},2));
for i = 1:size(m{1},2)
    % non-linear triangulationfor the i-th point
    M(:,i) = triang_nonlin(M0(:,i),P(vis(i,:)),num2cell(squeeze(points(:,i,vis(i,:))),1));
end

end


