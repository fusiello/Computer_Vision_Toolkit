function wireframe(vertices,edges,varargin)
    % draw a wireframe model with vertices and edges
    
    % Author: M. Farenzena, modified by A. Fusiello, 2007 (andrea.fusiello@univr,it)
    
    if ~ishold
        hold on
        wasNotHeld = 1;
    else
        wasNotHeld = 0;
    end
    
    x = vertices(:,1);
    y = vertices(:,2);
    if size(vertices,2) > 2
        z = vertices(:,3);
    end
    
    for i=1:size(edges,1)
        if size(vertices,2) > 2
            p1 = [x(edges(i,1)) y(edges(i,1)) z(edges(i,1))];
            p2 = [x(edges(i,2)) y(edges(i,2)) z(edges(i,2))];
            plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)], varargin{:});
        else
            p1 = [x(edges(i,1)) y(edges(i,1)) ];
            p2 = [x(edges(i,2)) y(edges(i,2)) ];
            plot([p1(1) p2(1)], [p1(2) p2(2)], varargin{:});
        end
    end
    
    if wasNotHeld
        hold off
    end
 
end

