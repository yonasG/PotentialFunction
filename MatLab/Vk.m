function [Vki] = Vk(k,i)
% Calculates the Legendre coefficients of the function f(t,i)
    f = @f1;
    Vki = (k+0.5).*integral(@(t) sin(t).*legendreP(k,cos(t)).*f(t,i),0,pi);
end 