function M  = triang_disp( u,v, d,  P1, P2 )
    %TRIANG_DISP Triangulation from disparity
    
    [K1,R1,t1] = krt(P1); [K2,R2,t2] = krt(P2);
    
    % baseline
    c1 = - R1'*t1; c2 = - R2'*t2;
    t = c1(1)-c2(1);
    
    T = [ 1/K1(1,1), 0,      0,         -K1(1,3)/K1(1,1)
        0,       1/K1(2,2),  0,         -K1(2,3)/K1(2,2)
        0,           0,      0,                 1
        0,           0,  1/(K1(1,1)*t), (K1(1,3)-K2(1,3))/(K1(1,1)*t)];
 
    
    M = htx(T, [u(:)';v(:)';d(:)']) ;
end







    % From symbolic with d=x2-x1 and t translation in camera 2 when
    % camera 1 is [I,0]
    %
    %
    % [ 1/au,    0,        0,            -u0/au]
    % [    0, 1/av,        0,            -v0/av]
    % [    0,    0,        0,                 1]   
    % [    0,    0, 1/(au*t), (u0 - u02)/(au*t)]

