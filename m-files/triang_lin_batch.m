
function M =triang_lin_batch(P, m, vis)
 %TRIANG_LIN_BATCH  Triangulation for n points in multiple images

% if no visibility, full visibility
if nargin < 4
    vis  = true(size(m{1},2),length(P));
end

% reshape the points cell matrix
    points  = cat(3,m{:});
    
M = NaN(3,size(m{1},2));
for i = 1:size(m{1},2)
    % triangulationfor the i-th point
      M(:,i) = triang_lin(P(vis(i,:)),num2cell(squeeze(points(:,i,vis(i,:))),1));
end

end


% function M = triang_lin_batch(P, m, vis)
%     %TRIANG_LIN_BATCH  Triangulation for n points in multiple images
%     
%     % if no visibility, full visibility
%      if nargin < 3
%         vis  = true(size(m{1},2),length(P));
%     end
%     
%     M=NaN(3,size(m{1},2));
%     for i = 1:size(m{1},2)
%         % triangulationfor the i-th point
%         L=[];
%         % stack equations for every camera
%         for j=1:length(P)
%             if vis(i,j)
%                 L = [L; skew([m{j}(:,i);1])*P{j} ]  ;
%             end
%         end
%         [~,~,V] = svd(L);
%         M(:,i)  = V(1:3,end)./V(4,end);
%     end
% end
