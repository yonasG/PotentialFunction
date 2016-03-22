function [out] = f3(t,i)
% Evaluates the function r^i for a given i   
    r = @distance;
    out =  r(t).^i;
end 