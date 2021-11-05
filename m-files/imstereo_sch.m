function [dmap,cost] = imstereo_sch(imL,imR,drange,ws)
    %IMSTEREO_SCH Stereo block-matching with SCH (Census)
    
    % compute Census transform 5x5
    fun = @(x) bin2dec(num2str(x(:)'>x(13)));
    imL = uint32(nlfilter(imL,[5 5],fun));
    imR = uint32(nlfilter(imR,[5 5],fun));
    
    dmin=drange(1); dmax=drange(2);
    sch=ones([size(imL),dmax-dmin+1])*Inf;
    for d = 0:dmax-dmin
        imR_d = circshift(imR,[0, -(dmin+d)]); % shift
        imR_d(:, end-(dmin+d):end) = 0;
        z = bitxor(imL, imR_d);
        ch = 0; % count the "1" in the xor
        for i=1:32 % 32 bits
            ch = ch + bitget(z,i);
        end
        sch(:,:,d+1)=filter2(ones(ws),ch);
    end
    [cost,dmap]=min(sch,[],3);
    dmap=dmap+dmin-1;
end

%% Warning: SCH needs 'nlfilter' from the Image Processing Toolbox