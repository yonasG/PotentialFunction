% PotentialFunction.m  
% Yonas Getachew
% 2/27/16
% This program calculates the coefficients of a real space potential that 
% corresponds to a Landau level mixing (LLM) pseudopotential in graphene
%
% Inputs: 
%       n = Landau level. 
%       Q = monopole strength
%       kappa = LLM parameter       
%       Vnm1 and dVnm1 = LLM and correction pseudopotentials
%       Type of correction function to add to the Coulomb potential
%    
% Outputs: 
%       coefficients for the expansion of the correection function
%       plot of the real space potential for varying kappa    
%*************************************************************************
clear all
clc
format long

% Define variables:
n = 0;              % Landau level
Q = 6.5;            % monopole strength
l = Q+n;            % angular monentum
kappa=1;

% Calls the function VCoul to caclculate the Coulomb pseudopotential V(L,Q)
Vc = zeros(2*l+1,1);
for L=1:2*l+1;
    Vc(L) = VCoul(L-1,Q);
end

fprintf('\n n = %.2f, Q = %.2f, kappa = %.2f, \n', n, Q, kappa)

% **********************************************************************
% Input: LLmxing pseudopotentials in n=1
Vnm1 = [0.829596; 0.599088; 0.469699; 0.372646; 0.326480; 0.297444;0.277387;...
        0.262864; 0.252102; 0.244088; 0.238205; 0.234064; 0.231421; 0.230133];

dVnm1 = [-0.4165; -0.1572; -0.0143; 0.0396; 0.0576; 0.0631; ... 
          0.0634; 0.0641; 0.0638; 0.0636; 0.0636; 0.0636; 0.0636; 0.0636];
      
% Constuct the matrix Vnm = Vnm1 + kappa*dVnm1 - V0
Vnm = zeros(2*l+1,1);
for L=1:2*l+1;
    Vnm(L) = Vnm1(L) + kappa.*dVnm1(L) - Vc(L);
end 
display(Vnm)

% **********************************************************************
% To check if the Legendre coefficients are being calculated correctly
% we find Vk of the Coulomb potential, which should be 1 for all k.
V0 = zeros(2*l+1,1);
for k=1:2*l+1;
    V0(k) = Vk0(k-1);
end

% **********************************************************************
% Crate the Vki(L,i) matrix
Vki = zeros(2*l+1,2*l+1);
for i=1:2*l+1;
    for k=1:2*l+1;
        Vki(k,i)= Vk(k-1,i);
        %fprintf('i, k, Vki(k,i) %.5f,  %.5f, %.5f \n', i, k, Vki(k,i))
    end
end
%display(Vki)

% **********************************************************************
% Crate the ALk(L,k) matrix using the Ak(Q,n,L,k) matrix
ALk = zeros(2*l+1,2*l+1);
for k=1:2*l+1;
    for L=1:2*l+1;
        ALk(L,k) = Ak(Q,n,L-1,k-1);
    end
end
%display(ALk)

% **********************************************************************
% multiply the matrices Vki and ALk to find dV
dV = mtimes(Vki,ALk);
%display(dV)
dV = dV - 0.001*eye(2*l+1);

% **********************************************************************
% Solve the system dV*C = Vnm for C using regular SVD method
[u, s, vT] = svd(dV);

% temp = u.'*Vnm;
% C = zeros(2*l+1,1);
% for k=1:2*l+1
%     C(k,:) = u(k).*temp(k).*1./s(i);
% end

C = vT*s*u.'*Vnm;

%C= vT*((u'*Vnm)./diag(s));

% check that A = U*S*VT
check1 = norm(dV - u*s*vT)/norm(dV);
% check the orthogonality of the right eigenvectors
check2 = norm(eye(2*l+1) - transpose(u)*u);
% check the orthogonality of the left eigenvectors
check3 = norm(eye(2*l+1) - vT*transpose(vT));
% check the difference between b and Ax
check4 = norm(Vnm - mtimes(dV,C))/norm(Vnm);

% **********************************************************************
fprintf('\n For kappa = %.2f The expansion coefficients are: \n', kappa)
display(C)

fprintf('\nThe three checks are: \n')
fprintf('Check 1 = %5.10E \n',check1)
fprintf('Check 2 = %5.10E \n',check2)
fprintf('Check 3 = %5.10E \n',check3)
fprintf('The value of |b-Ax|/|b| is: %5.10E \n', check4)


% C = [3.428269425656438;
%    2.230787865365076;
%    0.784111830735872;
%    0.136006548018487;
%    0.018996506538296;
%   -0.000827332978129;
%    0.001031331467866;
%    0.000049015670725;
%    0.000406417183675;
%    0.000395488035722;
%    0.000501757729597;
%    0.000644681334506;
%    0.000937595416990;
%    0.002125814184722];

%C = [1;1;1;1;1;1;1;1;1;1;1;1;1];
   
   
% **********************************************************************
% Test: this part checks if the pseudopotentials of the output potential
% function match the original pseudopotentials

% Finds the legendre coefficients of the new potential function
% for k=1:2*l+1
%      F=@potential;
%      Pcoeff = (k+0.5).*integral(@(t) sin(t).*legendreP(k,cos(t)).*F(t,C,l),0,pi);
%      %fprintf('\nThe %.2f th Legendre coefficeints is: %5.10f \n',k, Pcoeff)
% end

% Given the legendre coefficeints finds the pseudopotential VL
VL = zeros(2*l+1,1);
for j=1:2*l+1
    VLi = zeros(2*l+1,1);
    for k=1:2*l
        %fprintf('\nMade it in k loop, %.2f \n', k)
        VLi(k) = Vm(Q,n,j-1,k-1,C(k));
        %fprintf('VL(k) = %.15f \n', VLi(k))
    end
    VL(j) = sum(VLi); %/sqrt(Q);
end

display(VL)

% **********************************************************************
% plot the function V(r) = 1/r + c_i f_i(r) in the range 0<x<pi
r = @distance;
x = 0:.00314149:2*pi;
[p,q] = size(x);
fi= zeros(2*l+1,q);
for i=1:2*l+1
    for j=1:q
        fi(i,j) = C(i).*r(x(j)).^i.*exp(-r(x(j)));%   f1(x(j),i);
    end 
end

Coulomb = 1./r(x);
V = [Coulomb;fi];
Vr = sum(V);

% Plot the potential function
figure(1)
plot(x, Coulomb, 'b', x, Vr, 'r-.', 'LineWidth', 2)

% Set the axis limits
%axis([0 pi -10 30])

% Add title and axis labels
xlabel('Chord distance r')
ylabel('V(r)')

% Add a legend 
legend('Coulomb','Coulomb + f(r)')


