function [ pot ] = potential(t,C,l)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%     [n,m] = size(t);
%     fprintf('\n Size of t = %.2f by %.2f ', n, m)
%     fprintf('\nInside potential t = %.4f ', t)
%     fprintf('\ninside potential l = %.3f ', l)
    r = @distance;
    f = @f1;
    for i=1:2*l+1
        pot = 1./r(t) + C(i).*f(t,i);
    end
    
end

