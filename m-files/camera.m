function G = camera( cop, look, up )
%CAMERA Camera exterior orientation given COP,LOOK and UP
    
    R(3,:) = (look(:) - cop(:))';               % Z axis
    R(3,:) = R(3,:)/norm(R(3,:));
    R(1,:) = cross((up(:) - cop(:))',R(3,:))';  % X axis 
    R(1,:) = R(1,:)/norm(R(1,:));
    R(2,:) = cross(R(3,:),R(1,:));              % Y axis
    
    G = [R, -R*cop(:) ];
end



