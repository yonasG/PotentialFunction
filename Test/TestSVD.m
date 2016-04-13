% This program solve an example 3x3 system using the SVD method
clear all
clc

A =[3 2 -1; 2 -2 4; -1 0.5 -1];
b = [1;-2;0];

% with linsolve the solution comes out close to to 1, -2, -2
x = linsolve(A,b);
display(x)

% To solve this system using SVD we first decompose A into 
[u, s, vT] = svd(A);

display(u)
display(s)
display(vT)

[m,n] = size(A);
display(m)

y = vT*((u'*b)./diag(s));

% y =zeros(m,1);
% for i=1:m
%     %coeff = dot(vT(:,i),b)/s(i,i);
%     coeff = dot(u(:,i), vT(:,i))/s(i,i);
%     y = coeff*b;
%     display(coeff)
% end

display(y)