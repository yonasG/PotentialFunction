function [out] = Vm(Q,n,L,k, C)
% Calculates the pseudopotential of a function of r
% calls the functions Wigner3j(j1,j2,j3,m1,m2,m3) and Wigner6j(j1,j2,j3,J1,J2,J3)
    l = Q+n;
    %f = @potential;
    %out = f(C,t,l).*(-1)^(2*Q+L).*(2*l+1).^2.*Wigner6j(L,l,l,k,l,l).*Wigner3j(l,k,l,-Q,0,Q);
    out = C.*(-1)^(2*Q+L).*(2*l+1).^2.*Wigner6j(L,l,l,k,l,l).*Wigner3j(l,k,l,-Q,0,Q).^2;
    %out = (-1)^(2*Q+L).*(2*l+1).^2.*Wigner6j(L,l,l,k,l,l).*Wigner3j(l,k,l,-Q,0,Q).^2;
end

