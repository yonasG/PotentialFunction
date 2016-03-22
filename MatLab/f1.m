function [out] = f1(t,i)
% Evaluates the function r^i*e^-r for a given i
% To control from main write an if statement whicih picks a specific f
    r = @distance;

    out =  r(t).^i.*exp(-r(t));
    %out = exp(-r(t).^2).*r(t).^(2*i);
    %out =  sin(r(t)).*i;
    %out =  exp(-r(t)).*r(t).^i;
end 