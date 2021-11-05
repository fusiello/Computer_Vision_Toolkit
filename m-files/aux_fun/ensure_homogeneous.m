function m = ensure_homogeneous(m)
    %ENSURE_HOMOGENEOUS Make sure m is homogeneous 2D
   
    % do not work on 3D points
     if size(m,1) < 3
        m = [m;ones(1,size(m,2))];
     end  
end

