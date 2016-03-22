function [out] = f4(t,i)
% Evaluates the function r^i for a given i   
    r = @distance;
    out = cos(t).*r(t).^i;
end 