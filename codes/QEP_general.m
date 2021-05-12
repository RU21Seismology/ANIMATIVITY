function [Slowness, Polarization]=QEP_general(C,rho)
%% Documentation
% Function to solve the eigenvalue problem
% Added by Xiaoran Chen on 10/23/2020
% Matlab R2016a

global p 

R=[C(5,5),C(4,5),C(3,5);...
   C(4,5),C(4,4),C(3,4);...
   C(3,5),C(3,4),C(3,3)];

S=[2*C(1,5),C(1,4)+C(5,6),C(5,5)+C(3,1);...
   C(1,4)+C(5,6),2*C(4,6),C(4,5)+C(3,6);...
   C(5,5)+C(3,1),C(4,5)+C(3,6),2*C(3,5)];

T=[C(1,1),C(1,6),C(1,5);...
   C(1,6),C(6,6),C(6,5);...
   C(1,5),C(6,5),C(5,5)];

T_input=pinv(R)*(p*p.*T-rho.*diag([1 1 1]));
S_input=p*pinv(R)*S;
I_input=diag([1 1 1]);

%% Calculate polyeig
[eigvector, eigvalue]=polyeig(T_input,S_input,I_input);

[Slowness,Inx]=sort(eigvalue,'descend');
Polarization=eigvector;
for i=1:6
    Polarization(:,i)=eigvector(:,Inx(i));
end