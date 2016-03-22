function [out] = f2(t,i)
% Evaluates the function r^2i* e^(-r^2) for a given i   
    r = @distance;
    out = exp(-r(t).^2).*r(t).^(2*i);
end 