function [dmap,cost] = imstereo_ssd(imL,imR,drange,ws)
    %IMSTEREO_SSD Stereo block-matching with normalized SSD

    dmin=drange(1); dmax=drange(2);
    s1 = filter2(ones(ws),imL.^2);
    ssd = ones([size(imL),dmax-dmin+1])*Inf;
    for d = 0:dmax-dmin
        imR_d = circshift(imR,[0, -(dmin+d)]); % shift
        imR_d(:, end-(dmin+d):end) = 0;
        sd = (imL-imR_d).^2; % squared differences
        s2 = filter2(ones(ws),imR_d.^2);
        ssd(:,:,d+1) = filter2(ones(ws),sd)./sqrt(s1.*s2); % SSD
    end
    [cost,dmap]=min(ssd,[],3);
    dmap=dmap+dmin-1;
end
