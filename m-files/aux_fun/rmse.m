function r = rmse( x )
    %RMSE Root mean square error
     r = sqrt(sum(x(:).^2)/length(x)); 
end

