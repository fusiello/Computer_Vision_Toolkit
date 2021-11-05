function [v,J] = ang2cart( p )
    %ANG2CART Angular to cartesian coordinate conversion
    
    % give the cartesian coordinates of a point on the unit sphere
    % represented by means of its angular coordinates
    % supported dimensions: 2 (1 angle), 3 (2 angles) and 4 (3 angles)
    
    if length(p) == 1
        v = [cos(p(1));sin(p(1))];
        J = [-sin(p(1)); cos(p(1))];
    elseif length(p) == 2
        v = [cos(p(1));
            sin(p(1)) * cos(p(2));
            sin(p(1)) * sin(p(2))] ;
        J = [        -sin(p(1)),                0
            cos(p(1))*cos(p(2)), -sin(p(1))*sin(p(2))
            cos(p(1))*sin(p(2)),  cos(p(2))*sin(p(1))];
    elseif length(p) == 3
        v = [cos(p(1));
            sin(p(1)) * cos(p(2));
            sin(p(1)) * sin(p(2)) * cos(p(3));
            sin(p(1)) * sin(p(2)) * sin(p(3))] ;
        J = [                -sin(p(1)),                       0,                        0
            cos(p(1))*cos(p(2)),        -sin(p(1))*sin(p(2)),                        0
            cos(p(1))*cos(p(3))*sin(p(2)), cos(p(2))*cos(p(3))*sin(p(1)), -sin(p(1))*sin(p(2))*sin(p(3))
            cos(p(1))*sin(p(2))*sin(p(3)), cos(p(2))*sin(p(1))*sin(p(3)),  cos(p(3))*sin(p(1))*sin(p(2))];
    else
        error('ang2cart: dimension > 4 not supported')
    end
end



