function [AkL] = Ak(Q,n,L,k)
% Calculates the matrix element A(k,L).
% calls the functions Wigner3j(j1,j2,j3,m1,m2,m3) and Wigner6j(j1,j2,j3,J1,J2,J3)
    l = Q+n;
    AkL = (((-1)^(2*Q+L).*(2*l+1).^2).*Wigner6j(L,l,l,k,l,l).*Wigner3j(l,k,l,-Q,0,Q).^2)./sqrt(Q);
end 