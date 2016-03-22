function [VC] = VCoul(L,Q)
% Calculates the Coulumb pseudopotential
    VC = (2.*nchoosek(4*Q-2*L,2*Q-L).*nchoosek(4*Q+2*L+2,2*Q+L+1))/(nchoosek(4*Q+2,2*Q+1).^2);
end 