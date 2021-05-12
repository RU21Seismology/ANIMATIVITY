function [TU1, RUD1]=InfRecursive(E1, E2, w, z0, z1, z, D1,D2, Rud0)
%% Documentation
% Function to calculate the generalized transmission and reflection coefficients 
% Added by Xiaoran Chen on 08/16/2020
% Matlab R2016a

%% Edition starts from here
E1_11=E1(1:3,1:3);
E1_12=E1(1:3,4:6);
E1_21=E1(4:6,1:3);
E1_22=E1(4:6,4:6);

E2_11=E2(1:3,1:3);
E2_12=E2(1:3,4:6);
E2_21=E2(4:6,1:3);
E2_22=E2(4:6,4:6);

% Build up the reflection and transmission coefficients
% interface 1
delta1_d_z=[exp(-1i*w*D1(1)*(z-z0)) exp(-1i*w*D1(2)*(z-z0)) exp(-1i*w*D1(3)*(z-z0))];
delta2_u_z=[exp(-1i*w*D2(4)*(z-z1)) exp(-1i*w*D2(5)*(z-z1)) exp(-1i*w*D2(6)*(z-z1))];

deltaD_RT=diag([delta1_d_z delta2_u_z]);
RT=[E2_11 -E1_12;E2_21 -E1_22]\[E1_11 -E2_12;E1_21 -E2_22]*deltaD_RT;

Td1=RT(1:3,1:3);
Rud1=RT(1:3,4:6);
Rdu1=RT(4:6,1:3);
Tu1=RT(4:6,4:6);

RUD0=Rud0;
TU1=(diag([1 1 1])-Rdu1*RUD0)\Tu1;
RUD1=Td1*RUD0*TU1+Rud1;






