function [VkCoul] = Vk0(k)
% Calculates the Legendre coefficients of the Coulumb potential
    r = @distance;
    VkCoul = (k+0.5)*integral(@(t) sin(t).*legendreP(k,cos(t))./r(t),0,pi);
end 